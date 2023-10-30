
use crate::block_gen::{PhaseBlock, filter_out_alignment_record};
use crate::data_types::read_segments::ReadSegment;
use crate::data_types::reference_genome::ReferenceGenome;
use crate::data_types::variants::{Variant, VariantType};
use crate::wfa_graph::{NodeAlleleMap, WFAGraph, WFAResult};
use crate::writers::phase_stats::ReadStats;

use bio::data_structures::interval_tree::IntervalTree;
use log::{debug, trace, warn};
use rust_htslib::bam;
use rustc_hash::FxHashMap as HashMap;
use simple_error::bail;
use std::path::{Path, PathBuf};

/// Loads up all the reads in a particular phase region and converts them into their variant representation.
/// This version uses local re-alignment to parse the alleles.
/// Returns an interval tree containing all reads to use for phasing, a second tree containing extra reads that *can* be phased but didn't match our criteria,
/// and statistics from loading the reads.
/// # Arguments
/// * `phase_problem` - the phase block we are loading data for
/// * `bam_paths` - the BAM files to parse, must be indexed
/// * `reference_filename` - the reference fasta file
/// * `variant_calls` - the variants used to convert full reads into haplotype observations (`ReadSegment`)
/// * `min_matched_alleles` - the minimum number of identified alleles required for a read to be included 
/// * `min_mapq` - the minimum MAPQ to consider a read
#[allow(clippy::type_complexity)]
pub fn load_read_segments(
    phase_problem: &PhaseBlock, bam_paths: &[PathBuf], reference_filename: &Path,
    variant_calls: &[Variant], min_matched_alleles: usize, min_mapq: u8
) -> Result<(IntervalTree<usize, ReadSegment>, IntervalTree<usize, ReadSegment>, ReadStats), Box<dyn std::error::Error>> {
    use rust_htslib::bam::Read;
    use rust_htslib::bam::ext::BamRecordExtensions;
    
    let num_variants: usize = variant_calls.len();
    let mut read_groups: HashMap<String, Vec<ReadSegment>> = Default::default();

    // stats we track
    let mut num_reads: u64 = 0;
    let mut skipped_reads: u64 = 0;
    let mut num_alleles: u64 = 0;
    let mut exact_matches: [u64; VariantType::Unknown as usize + 1] = [0; VariantType::Unknown as usize + 1];
    let mut inexact_matches: [u64; VariantType::Unknown as usize + 1] = [0; VariantType::Unknown as usize + 1];
    let mut failed_matches: [u64; VariantType::Unknown as usize + 1] = [0; VariantType::Unknown as usize + 1];
    let mut allele0_matches: [u64; VariantType::Unknown as usize + 1] = [0; VariantType::Unknown as usize + 1];
    let mut allele1_matches: [u64; VariantType::Unknown as usize + 1] = [0; VariantType::Unknown as usize + 1];
    
    for bam_filename in bam_paths.iter() {
        let mut bam_reader = bam::IndexedReader::from_path(bam_filename)?;
        bam_reader.set_reference(reference_filename)?;
        bam_reader.fetch((phase_problem.get_chrom(), phase_problem.get_start(), phase_problem.get_end()+1))?;
        
        for read_entry in bam_reader.records() {
            let mut read = read_entry?;
            
            //make sure we care about the alignment
            if filter_out_alignment_record(&read, min_mapq) {
                continue;
            }

            //build out the cigar info
            read.cache_cigar();
            
            //build a lookup from reference coordinate -> sequence coordinate
            let mut coordinate_lookup: HashMap<i64, i64> = Default::default();
            let min_position = read.pos();
            let mut max_position = read.pos();
            for bp in read.aligned_pairs() {
                let segment_index = bp[0];
                let ref_index = bp[1];
                coordinate_lookup.insert(ref_index, segment_index);
                max_position = max_position.max(ref_index);
            }
            assert!(max_position >= min_position);

            // max_position is the last one that we found, so add +1 to include in the range
            let aligned_range = min_position..(max_position+1);

            //.seq() returns Seq<'_> type, but we should just full decode
            let read_sequence: Vec<u8> = read.seq().as_bytes();
            let read_qualities: &[u8] = read.qual();
            assert_eq!(read_sequence.len(), read_qualities.len());

            //we will populate these with the variant level info
            let mut alleles: Vec<u8> = Vec::<u8>::with_capacity(num_variants);
            let mut quals: Vec<u8> = Vec::<u8>::with_capacity(num_variants);
            let mut num_overlaps: usize = 0;
            let mut last_deletion_end: usize = 0;

            for variant in variant_calls.iter() {
                /*
                - We need to split on small variants and SVs
                - for small variants, do what we normally do; it may be worth seeing if the method we create for SVs will help with this other mode fails though
                - for SVs, check if the read FULLY spans the locus; if so, check how much sequence is inserted/deleted in the region and turn that into an allele
                - TODO: for SVs, what if it doesn't, can we use clipping somehow?
                */

                trace!("{:?}", variant);
                let variant_pos: i64 = variant.position();
                let variant_type: VariantType = variant.get_type();
                let vt_index = variant_type as usize;

                // regardless of variant type, we MUST populate these in the following branching logic
                let mut allele: u8;
                let qual: u8;
                let exact_allele: bool;
                let overlaps_allele: bool;

                if variant.is_ignored() {
                    // this variant is one marked to ignored, lets set it to undefined as opposed to ambiguous
                    trace!("\tMarking as undefined allele because it is flagged to be ignored");
                    allele = 3;
                    qual = 0;
                    exact_allele = false;
                    overlaps_allele = false;
                } else if variant_pos < last_deletion_end as i64 {
                    // check if this is within a region we have decided is a deleted
                    trace!("\tMarking as unknown allele because it overlaps detected SV deletion");
                    // if the 0-allele is reference, mark as 0, else mark as ambiguous because it's multi-allelic call
                    allele = 2;
                    qual = 0;
                    exact_allele = false;
                    overlaps_allele = true;
                } else {
                    match variant_type {
                        VariantType::Snv | 
                        VariantType::Insertion | 
                        VariantType::Deletion |
                        VariantType::Indel |
                        VariantType::SvInsertion |
                        VariantType::TandemRepeat => {
                            // we need these to build coordinate ranges
                            let ref_allele_len: usize = variant.get_ref_len();
                            let prefix_len: usize = variant.get_prefix_len();
                            let postfix_len: usize = variant.get_postfix_len();
                            
                            // coordinate ranges we care about
                            let first_start_coordinate: usize = variant_pos as usize - prefix_len;
                            let last_start_coordinate: usize = variant_pos as usize + 1; // add one because we want to include variant_pos
                            let first_end_coordinate: usize = variant_pos as usize + ref_allele_len;
                            let last_end_coordinate: usize = variant_pos as usize + ref_allele_len + postfix_len + 1; // add one for same reason as above

                            // first, try to find the closest start
                            let mut opt_closest_start: Option<usize> = None;
                            for sc in (first_start_coordinate..last_start_coordinate).rev() {
                                if let Some(&si) = coordinate_lookup.get(&(sc as i64)) {
                                    opt_closest_start = Some(si as usize);
                                    break;
                                }
                            }

                            // now the closest end
                            let mut opt_closest_end: Option<usize> = None;
                            for ec in first_end_coordinate..last_end_coordinate {
                                if let Some(&ei) = coordinate_lookup.get(&(ec as i64)) {
                                    opt_closest_end = Some(ei as usize);
                                    break;
                                }
                            }

                            // now find the best start coordinate with constraints
                            let mut start_coordinate: Option<usize> = None;
                            let mut start_clip: usize = 0;
                            let mut end_coordinate: Option<usize> = None;
                            let mut end_clip: usize = 0;

                            if let (Some(closest_start), Some(closest_end)) = (opt_closest_start, opt_closest_end) {
                                for sc in first_start_coordinate..last_start_coordinate {
                                    // always increment this
                                    start_clip += 1;

                                    if let Some(&segment_index) = coordinate_lookup.get(&(sc as i64)) {
                                        // check if it's too far away
                                        if closest_start - segment_index as usize > 2*prefix_len {
                                            continue;
                                        }

                                        // we found a start coordinate
                                        start_coordinate = Some(segment_index as usize);

                                        // now try to find an end coordinate also
                                        for ec in (first_end_coordinate..last_end_coordinate).rev() {
                                            // always increment this
                                            end_clip += 1;

                                            if let Some(&next_index) = coordinate_lookup.get(&(ec as i64)) {
                                                // check if it's too far away
                                                if next_index as usize - closest_end > 2*postfix_len {
                                                    continue;
                                                }
                                                
                                                // we found an end coordinate also
                                                end_coordinate = Some(next_index as usize);
                                                break;
                                            }
                                        }
                                        break;
                                    }
                                }
                            } else {
                                // the closest ones failed, we won't succeed here either
                            }

                            if let Some(ss) = start_coordinate {
                                if let Some(se) = end_coordinate {
                                    trace!("\t{}..{} = {:?} {:?}; next = {} {}", ss, se, &read_sequence[ss..se], &read_qualities[ss..se], read_sequence[se], read_qualities[se]);
                                    
                                    let edit_distance: usize;
                                    allele = variant.match_allele(&read_sequence[ss..se]);
                                    if allele == 2 {
                                        // no exact match, do inexact matching
                                        (allele, edit_distance, _) = variant.closest_allele_clip(&read_sequence[ss..se], start_clip - 1, end_clip - 1);
                                        exact_allele = false;
                                    } else {
                                        edit_distance = 0;
                                        exact_allele = true;
                                    }

                                    // this approach uses harmonic mean of base quality
                                    // * no ED penalty
                                    // * weighted - the same weight factors are applied, but in a down-weighting approach (this is because we don't want to exceed u8::MAX)
                                    let divisor_multiplier: i64 = match variant_type {
                                        // these are down-weights for local mode
                                        // SNV has highest confidence here
                                        VariantType::Snv => 1,

                                        // indels tend to be pretty bad
                                        VariantType::Deletion |
                                        VariantType::Insertion |
                                        VariantType::Indel => 4,

                                        // SVs are generally worst of all in local mode
                                        VariantType::SvDeletion |
                                        VariantType::SvInsertion => 4,

                                        // we want tandem repeats to have a higher weight than our generic indels
                                        VariantType::TandemRepeat => 2,

                                        _ => {
                                            panic!("No implementation for matching {variant_type:?}");
                                        }
                                    };
                                    
                                    qual = ((
                                        (se - ss) as f64 /
                                        read_qualities[ss..se].iter()
                                            .map(|&q| 1.0f64 / q as f64)
                                            .sum::<f64>()
                                    ).max(4.0) / divisor_multiplier as f64) as u8;
                                    
                                    overlaps_allele = true;
                                    trace!("\tallele = {}, qual = {}, ED = {}", allele, qual, edit_distance);
                                } else {
                                    trace!("\tfailed allele match for ref extension");
                                    allele = 2;
                                    qual = 0;
                                    exact_allele = false;
                                    overlaps_allele = true;
                                }
                            } else {
                                //no overlap
                                if aligned_range.contains(&variant_pos) {
                                    trace!("\tOverlap, no position");
                                    overlaps_allele = true;
                                    allele = 2;
                                } else {
                                    // there is no alignment overlap
                                    trace!("\tNo overlap");
                                    overlaps_allele = false;
                                    allele = 3;
                                }
                                qual = 0;
                                exact_allele = false;
                            }
                        },
                        VariantType::SvDeletion => {
                            // we need these to build coordinate ranges
                            let ref_allele_len: usize = variant.get_ref_len();
                            
                            if aligned_range.contains(&variant_pos) {
                                // coordinate ranges we care about
                                let last_start_coordinate: usize = variant_pos as usize + 1; // add one because we want to include variant_pos
                                let first_end_coordinate: usize = variant_pos as usize + ref_allele_len;
                                if aligned_range.contains(&(first_end_coordinate as i64)) {
                                    // calculate how many bases we expect to see deleted
                                    let expected_deleted: usize = first_end_coordinate - last_start_coordinate;
                                    
                                    // now we need to move up and down until we find an anchor point
                                    let mut start_anchor: usize = last_start_coordinate;
                                    while !coordinate_lookup.contains_key(&(start_anchor as i64)) {
                                        if start_anchor <= aligned_range.start as usize {
                                            // fixes weird CIGARs where a mapping starts with non-matching types, e.g.:
                                            //   [SoftClip(3139), Del(798), Equal(4), ...
                                            warn!("Reached start of read ({}) without finding start_anchor, using POS ({}) instead.", std::str::from_utf8(read.qname()).unwrap_or("utf8-decode-error"), start_anchor);
                                            break;
                                        }
                                        start_anchor -= 1;
                                    }
                                    let mut end_anchor: usize = first_end_coordinate;
                                    while !coordinate_lookup.contains_key(&(end_anchor as i64)) {
                                        end_anchor += 1;
                                        if end_anchor >= aligned_range.end as usize {
                                            // we have not observed it, but this is a symmetrical handling of the weird CIGARs for the end
                                            warn!("Reached end of read ({}) without finding end_anchor, using max ({}) found instead.", std::str::from_utf8(read.qname()).unwrap_or("utf8-decode-error"), end_anchor);
                                            break;
                                        }
                                    }
                                    
                                    // count up the number of missing (i.e. deleted) based in the reference
                                    let mut deleted_count: usize = 0;
                                    for dc in start_anchor..end_anchor {
                                        if !coordinate_lookup.contains_key(&(dc as i64)) {
                                            deleted_count += 1;
                                        }
                                    }
                                    
                                    // it's possible to have more deleted bases than expected
                                    // assert!(expected_deleted >= deleted_count);
                                    
                                    // the quality if we have exactly the right number of deleted (or not deleted) bases
                                    let exact_allele_qual: f64 = 40.0;
                                    // divisor for downweighting quality values, this is higher than indels currently
                                    let deletion_factor: f64 = 8.0;
                                    // fixes the ratios that match REF or ALT here to: REF = [0, match_window_size); ALT = (1.0 - match_window_size, 1.0 + match_window_size)
                                    let match_window_size: f64 = 0.33;
                                    
                                    let deleted_ratio: f64 = deleted_count as f64 / expected_deleted as f64;
                                    if deleted_ratio < match_window_size {
                                        // mostly not deleted
                                        allele = 0;
                                        if deleted_ratio == 0.0 {
                                            // this is pretty unlikely
                                            qual = (exact_allele_qual / deletion_factor) as u8;
                                            exact_allele = true;
                                        } else {
                                            qual = ((-10.0 * deleted_ratio.log10()) / deletion_factor).max(1.0) as u8;
                                            exact_allele = false;
                                        }
                                    } else if (1.0 - deleted_ratio).abs() < match_window_size  {
                                        // mostly deleted and not over-deleted
                                        allele = 1;
                                        if deleted_ratio == 1.0 {
                                            // this is pretty unlikely
                                            qual = (exact_allele_qual / deletion_factor) as u8;
                                            exact_allele = true;
                                        } else {
                                            qual = ((-10.0 * (1.0 - deleted_ratio).abs().log10()) / deletion_factor).max(1.0) as u8;
                                            exact_allele = false;
                                        }

                                        // this is getting labeled a deletion, force anything overlapping it to be reference (because it isn't there)
                                        last_deletion_end = first_end_coordinate;
                                    } else {
                                        // ambiguous either because it's in between or over-deleted
                                        allele = 2;
                                        qual = 0;
                                        exact_allele = false;
                                    }
                                    overlaps_allele = true;
                                } else {
                                    // we have a partial overlap, but don't reach the far end
                                    // mirror what we do above by marking overlap as true but otherwise a failure to match
                                    allele = 2;
                                    qual = 0;
                                    exact_allele = false;
                                    overlaps_allele = true;
                                }
                            } else {
                                // we don't overlap the start
                                allele = 3;
                                qual = 0;
                                exact_allele = false;
                                overlaps_allele = false;
                            }
                        },
                        _ => {
                            panic!("Unhandled variant type: {variant_type:?}");
                        }
                    };
                }

                // gather stats on the match
                if overlaps_allele {
                    assert!(allele <= 2);
                    if allele == 2 {
                        failed_matches[vt_index] += 1;
                    } else {
                        if exact_allele {
                            exact_matches[vt_index] += 1;
                        } else {
                            inexact_matches[vt_index] += 1;
                        }
                        if allele == 0 {
                            allele0_matches[vt_index] += 1;
                        } else {
                            allele1_matches[vt_index] += 1;
                        }
                        num_overlaps += 1;
                        num_alleles += 1;
                    }
                } else {
                    assert_eq!(allele, 3);
                }

                // no matter what, we push these now
                alleles.push(allele);
                // make sure the quality is always at least 1
                quals.push(qual.max(1));
            }
            assert_eq!(num_variants, alleles.len());
            assert_eq!(num_variants, quals.len());
            trace!("All alleles {:?}\n", alleles);

            if num_overlaps > 0 {
                let read_name: String = String::from_utf8(read.qname().to_vec()).unwrap();
                let read_group: &mut Vec<ReadSegment> = read_groups.entry(read_name.clone()).or_insert(vec![]);
                read_group.push(ReadSegment::new(read_name, alleles, quals));
            } else {
                // this one has no overlaps, so it's just a skipped read
                skipped_reads += 1;
            }
        }
    }

    // now collapse all the reads, but only keeping those with at least 2 things set
    let mut read_segments: IntervalTree<usize, ReadSegment> = IntervalTree::new();
    let mut phasable_segments: IntervalTree<usize, ReadSegment> = IntervalTree::new();
    for (_qname, read_group) in read_groups.iter() {
        let collapsed_read: ReadSegment = ReadSegment::collapse(read_group);
        let num_set: usize = collapsed_read.get_num_set();
        if num_set >= min_matched_alleles {
            let segment_range = collapsed_read.get_range();
            read_segments.insert(segment_range, collapsed_read);
            num_reads += read_group.len() as u64;
        } else {
            skipped_reads += read_group.len() as u64;
            if num_set > 0 {
                // even though this won't be used for phasing, it CAN be phased
                let segment_range = collapsed_read.get_range();
                phasable_segments.insert(segment_range, collapsed_read);
            }
        }
    }
    
    // sanity check this; this was before we started making sure failed alleles only applied if the mapping overlapped
    // assert_eq!(num_alleles, (num_reads + skipped_reads) * (num_variants as u64));
    let segment_stats = ReadStats::new(
        num_reads, skipped_reads, num_alleles, 
        exact_matches, inexact_matches, failed_matches,
        allele0_matches, allele1_matches,
        false
    );
    debug!("Read segment stats: {:?}", segment_stats);
    
    Ok((read_segments, phasable_segments, segment_stats))
}

/// Loads up all the reads in a particular phase region and converts them into their variant representation.
/// This version uses global re-alignment to parse the alleles.
/// Returns an interval tree containing all reads to use for phasing, a second tree containing extra reads that *can* be phased but didn't match our criteria,
/// and statistics from loading the reads.
/// # Arguments
/// * `phase_problem` - the phase block we are loading data for
/// * `bam_paths` - the BAM files to parse, must be indexed
/// * `variant_calls` - the variants used to convert full reads into haplotype observations (`ReadSegment`)
/// * `hom_calls` - any homozygous variants within the region, these don't get phased but are useful for global realignment
/// * `reference_genome` - the reference genome sequences, required for this approach
/// * `min_matched_alleles` - the minimum number of identified alleles required for a read to be included 
/// * `min_mapq` - the minimum MAPQ to consider a read
/// * `max_runtime` - controls the allowed runtime of the global realignment
/// * `wfa_prune_distance` - maximum allowed distance a wavefront can lag; make smaller to reduce run-time at the cost of accuracy
#[allow(clippy::too_many_arguments)]
#[allow(clippy::type_complexity)]
pub fn load_full_read_segments(
    phase_problem: &PhaseBlock, bam_paths: &[PathBuf], variant_calls: &[Variant], hom_calls: &[Variant],
    reference_genome: &ReferenceGenome, min_matched_alleles: usize, min_mapq: u8,
    max_runtime: f32, wfa_prune_distance: usize
) -> Result<(IntervalTree<usize, ReadSegment>, IntervalTree<usize, ReadSegment>, ReadStats), Box<dyn std::error::Error>> {
    use rust_htslib::bam::Read;
    use rust_htslib::bam::ext::BamRecordExtensions;
    
    let chromosome: &str = phase_problem.get_chrom();
    let num_variants: usize = variant_calls.len();
    let mut read_groups: HashMap<String, Vec<ReadSegment>> = Default::default();

    // stats we track
    let mut num_reads: u64 = 0;
    let mut skipped_reads: u64 = 0;
    let mut num_alleles: u64 = 0;
    let mut exact_matches: [u64; VariantType::Unknown as usize + 1] = [0; VariantType::Unknown as usize + 1];
    let mut inexact_matches: [u64; VariantType::Unknown as usize + 1] = [0; VariantType::Unknown as usize + 1];
    let mut failed_matches: [u64; VariantType::Unknown as usize + 1] = [0; VariantType::Unknown as usize + 1];
    let mut allele0_matches: [u64; VariantType::Unknown as usize + 1] = [0; VariantType::Unknown as usize + 1];
    let mut allele1_matches: [u64; VariantType::Unknown as usize + 1] = [0; VariantType::Unknown as usize + 1];
    
    let mut edit_distances: Vec<usize> = vec![];

    let block_cpu_time = cpu_time::ThreadTime::now();

    for bam_filename in bam_paths.iter() {
        let mut bam_reader = bam::IndexedReader::from_path(bam_filename)?;
        bam_reader.set_reference(reference_genome.filename())?;
        bam_reader.fetch((chromosome, phase_problem.get_start(), phase_problem.get_end()+1))?;
        
        for read_entry in bam_reader.records() {
            let time_elapsed: f32 = block_cpu_time.elapsed().as_secs_f32();
            if time_elapsed > max_runtime {
                bail!("max_runtime reached");
            }

            let mut read = read_entry?;
            
            //make sure we care about the alignment
            if filter_out_alignment_record(&read, min_mapq) {
                continue;
            }

            //build out the cigar info
            read.cache_cigar();
            
            //build a lookup from reference coordinate -> sequence coordinate
            let mut coordinate_lookup: HashMap<i64, i64> = Default::default();
            let mut min_position: i64 = i64::MAX;
            let mut max_position: i64 = i64::MIN;
            for bp in read.aligned_pairs() {
                let segment_index = bp[0];
                let ref_index = bp[1];
                coordinate_lookup.insert(ref_index, segment_index);
                min_position = min_position.min(ref_index);
                max_position = max_position.max(ref_index);
            }
            assert!(max_position >= min_position);

            // max_position is the last one that we found, so add +1 to include in the range
            let aligned_range = min_position..(max_position+1);

            //we will populate these with the variant level info
            let mut num_overlaps: usize = 0;
            let mut first_overlap: Option<usize> = None;
            let mut last_overlap: usize = 0;
            for (i, variant) in variant_calls.iter().enumerate() {
                let variant_pos: i64 = variant.position();
                if aligned_range.contains(&variant_pos) {
                    if first_overlap.is_none() {
                        first_overlap = Some(i);
                    }
                    last_overlap = i+1;
                    num_overlaps += 1;
                }
            }

            // if this mapping overlaps no alleles, then there's no reason to look at it anymore
            if num_overlaps == 0 {
                skipped_reads += 1;
                continue;
            }

            // convert into a non-option
            let first_overlap: usize = first_overlap.unwrap();
            assert_eq!(num_overlaps, last_overlap - first_overlap);

            // check for homozygous variants also
            let mut first_hom_overlap: Option<usize> = None;
            let mut last_hom_overlap: usize = 0;
            for (i, variant) in hom_calls.iter().enumerate() {
                let variant_pos: i64 = variant.position();
                if aligned_range.contains(&variant_pos) {
                    if first_hom_overlap.is_none() {
                        first_hom_overlap = Some(i);
                    }
                    last_hom_overlap = i+1;
                }
            }
            let first_hom_overlap: usize = first_hom_overlap.unwrap_or(0);

            // .seq() returns Seq<'_> type, but we should just full decode
            let read_sequence: Vec<u8> = read.seq().as_bytes();
            let read_qualities: &[u8] = read.qual();
            assert_eq!(read_sequence.len(), read_qualities.len());
            
            // these should always exist based on how we set it up
            let read_start: usize = *coordinate_lookup.get(&min_position).unwrap() as usize;
            let read_end: usize = *coordinate_lookup.get(&max_position).unwrap() as usize;

            // pull out the part of the read we're aligning against
            let read_align: &[u8] = &read_sequence[read_start..(read_end+1)];

            /*
            Current state:
            - we have the reference genome
            - we have the part of the read that aligns in `read_align`, the full read sequence in `read_sequence`
            - we have the indices of the first and last variant overlaps in `first_overlap` and `last_overlap`
            
            We need to populate:
            - alleles
            - quals
            - read stats (see below)

            Game plan:
            - construct a graph representing just this reference location + relevant alleles
            - while constructing, assign alleles to each new branch (it may be reference allele)
            -- IF you have multiple alleles starting at the same coordinate (e.g. identical call), then do not create an in-between node; this should resolve in the tie-breaking as "identical"
            -- so each branch should get a variant index + an allele assignment (0/1); reference alleles may end up with multiple 0 alleles in the event of multi-start
            - align the read via POA
            - look at the traversed nodes and copy the allele assignments; if anything is unassigned at the end, it gets 2; any with conflicting assignments get 2 also
            - update stats according to the assignments, we can't really do exact right now (maybe we can look at score deltas from one node to the next?)
            */
            let chrom_seq: &[u8] = reference_genome.get_full_chromosome(chromosome);
            
            // we need to also provide any preset alleles
            let start_time = std::time::Instant::now();
            let (wfa_graph, node_to_alleles): (WFAGraph, NodeAlleleMap) = 
                WFAGraph::from_reference_variants_with_hom(
                    chrom_seq, 
                    &variant_calls[first_overlap..last_overlap], // these are both range style indices
                    &hom_calls[first_hom_overlap..last_hom_overlap],
                    min_position as usize, 
                    max_position as usize + 1
                ).unwrap();
            
            let wfa_result: WFAResult = match wfa_graph.edit_distance_with_pruning(read_align, wfa_prune_distance) {
                Ok(wr) => wr,
                Err(e) => {
                    bail!(
                        "Encountered WFA error for mapping \"{}\" ({}:{}): {}", 
                        std::str::from_utf8(read.qname()).unwrap_or("QNAME_UTF8_ERROR"), 
                        chromosome, read.pos(), e
                    );
                }
            };
            
            debug!(
                "B#{} WFAGraph result ({}) => num_nodes: {}, read_len: {}, variant_overlaps: {}, edit_distance: {}", 
                phase_problem.get_block_index(), start_time.elapsed().as_secs_f32(), wfa_graph.get_num_nodes(), max_position-min_position+1, num_overlaps, wfa_result.score()
            );
            
            edit_distances.push(wfa_result.score());

            //we will populate these with the variant level info
            let mut alleles: Vec<u8> = vec![3; num_variants];
            for traversed_index in wfa_result.traversed_nodes().iter() {
                for &(var_index, allele_assignment) in node_to_alleles.get(traversed_index).unwrap_or(&vec![]).iter() {
                    let correct_index: usize = first_overlap+var_index;
                    if alleles[correct_index] == 3 {
                        alleles[correct_index] = allele_assignment;
                    } else if alleles[correct_index] != allele_assignment {
                        alleles[correct_index] = 2;
                    }
                }
            }
            
            // go through the result counting assigned and setting qualities
            let mut quals: Vec<u8> = vec![0; num_variants];
            for (i, a) in alleles.iter_mut().enumerate() {
                let variant_type: VariantType = variant_calls[i].get_type();
                let vt_index: usize = variant_type as usize;
                if *a == 3 {
                    // no overlaps for this allele
                } else if *a == 2 {
                    // overlaps, but ambiguous matching
                    failed_matches[vt_index] += 1;
                } else {
                    // we got a match, figure out the quality for it
                    quals[i] = match variant_type {
                        // these weights are up-weights for global re-alignments
                        // SNVs tend to always be the cleanest
                        VariantType::Snv => 8,

                        // these are probably the noisiest of the bunch
                        VariantType::Deletion |
                        VariantType::Insertion |
                        VariantType::Indel => 1,
                        
                        // these should be pretty high confidence because they have a lot of bases to make them work
                        VariantType::SvDeletion |
                        VariantType::SvInsertion => 2,
                        
                        // we want tandem repeats to have higher confidence than random indels
                        VariantType::TandemRepeat => 4,

                        _ => {
                            panic!("No implementation for matching {variant_type:?}");
                        }
                    };

                    // gather stats on the match
                    let exact_allele = false; // TODO: figure this out
                    if exact_allele {
                        exact_matches[vt_index] += 1;
                    } else {
                        inexact_matches[vt_index] += 1;
                    }
                    if *a == 0 {
                        allele0_matches[vt_index] += 1;
                    } else {
                        allele1_matches[vt_index] += 1;
                    }
                    num_alleles += 1;
                }
            } 

            // need to check what these were before
            assert_eq!(num_variants, alleles.len());
            assert_eq!(num_variants, quals.len());
            trace!("All alleles {:?}\n", alleles);

            let read_name: String = String::from_utf8(read.qname().to_vec()).unwrap();
            let read_group: &mut Vec<ReadSegment> = read_groups.entry(read_name.clone()).or_insert(vec![]);
            read_group.push(ReadSegment::new(read_name, alleles, quals));
        }
    }

    // now collapse all the reads, but only keeping those with at least 2 things set
    let mut read_segments: IntervalTree<usize, ReadSegment> = IntervalTree::new();
    let mut phasable_segments: IntervalTree<usize, ReadSegment> = IntervalTree::new();
    for (_qname, read_group) in read_groups.iter() {
        let collapsed_read: ReadSegment = ReadSegment::collapse(read_group);
        let num_set: usize = collapsed_read.get_num_set();
        if num_set >= min_matched_alleles {
            let segment_range = collapsed_read.get_range();
            read_segments.insert(segment_range, collapsed_read);
            num_reads += read_group.len() as u64;
        } else {
            skipped_reads += read_group.len() as u64;
            if num_set > 0 {
                // even though this won't be used for phasing, it CAN be phased
                let segment_range = collapsed_read.get_range();
                phasable_segments.insert(segment_range, collapsed_read);
            }
        }
    }
    
    // sanity check this; this was before we started making sure failed alleles only applied if the mapping overlapped
    // assert_eq!(num_alleles, (num_reads + skipped_reads) * (num_variants as u64));
    let segment_stats = ReadStats::new(
        num_reads, skipped_reads, num_alleles, 
        exact_matches, inexact_matches, failed_matches,
        allele0_matches, allele1_matches,
        true
    );
    debug!("Read segment stats: {:?}", segment_stats);
    debug!("Edit distances: {:?}", edit_distances);
    
    Ok((read_segments, phasable_segments, segment_stats))
}