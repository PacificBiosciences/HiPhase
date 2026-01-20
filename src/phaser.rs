
use crate::astar_phaser::{self, AstarConfig, AstarResult};
use crate::block_gen::{PhaseBlock, is_phasable_variant, get_variant_type};
use crate::data_types::read_segments::{AlleleType, ReadSegment};
use crate::data_types::reference_genome::ReferenceGenome;
use crate::data_types::variants::{IgnoredVariantReason, Variant, VariantType};
use crate::read_parsing::{self, GlobalRealignmentConfig};
use crate::variant_traversal::optimize_variant_order;
use crate::writers::phase_stats::{PhaseStats, ReadStats};

use bio::data_structures::interval_tree::IntervalTree;
use log::{debug, trace};
use priority_queue::PriorityQueue;
use rust_htslib::bcf;
use rust_htslib::bcf::record::GenotypeAllele;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use simple_error::{SimpleError, bail};
use std::cmp::{Ordering, Reverse};
use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;

/// Configuration for the phasing algorithm.
pub struct PhasingConfig {
    // variable related to read parsing
    /// The reference buffer
    pub reference_buffer: usize,
    /// The minimum number of matched alleles required to include a read
    pub min_matched_alleles: usize,
    /// The minimum MAPQ to include a read
    pub min_mapq: u8,
    /// The global realignment configuration
    pub global_realignment_config: Option<read_parsing::GlobalRealignmentConfig>,

    // variable related to the phasing algorithm
    /// if true, then we will reorder the variants before phasing
    pub optimize_variant_order: bool,
    /// The A* configuration
    pub astar_config: AstarConfig,

    // variables related to the post-phasing block splitting
    /// if true, then we will use the old linear block algorithm
    pub use_linear_block_algorithm: bool,
    /// the minimum number of connecting reads to use for the non-linearblock algorithm
    pub min_connecting_reads: u64,
}

/// Core function for loading variant calls from our VCF file and converting them into a `Variant` type.
/// # Arguments
/// * `region` - the phase block we need to load
/// * `vcf_paths` - the VCF files to load, must be zipped and indexed
/// * `reference_genome` - optional, the reference genome 
/// * `reference_buffer` - the number of nearby bases to try to use for local realignment
/// * `is_hom_allowed` - if true, then non-reference homozygous variants will also be loaded into the second return Vec
fn load_variant_calls(
    region: &PhaseBlock, 
    vcf_paths: &[PathBuf],
    reference_genome: &ReferenceGenome, 
    reference_buffer: usize, 
    is_hom_allowed: bool
) -> Result<(Vec<Variant>, Vec<Variant>), Box<dyn std::error::Error>> {
    use rust_htslib::bcf::Read;
    
    // short circuit because otherwise bcf can throw errors
    if region.get_num_variants() == 0 {
        return Ok((vec![], vec![]));
    }

    // initalize the variant queue with one variant from each VCF, any ties in position are broken by VCF input order
    let mut variant_queue: PriorityQueue<usize, (Reverse<i64>, Reverse<usize>)> = PriorityQueue::new();
    let mut vcf_readers: Vec<bcf::IndexedReader> = vcf_paths.iter()
        .map(|filename| bcf::IndexedReader::from_path(filename).unwrap())
        .collect();
    let mut vcf_iterators: Vec<_> = vec![];
    let mut sample_indices: Vec<usize> = vec![];
    let sample_name: &str = region.sample_name();

    // fetch the region for each VCF (if it exists)
    for (vcf_index, vcf_reader) in vcf_readers.iter_mut().enumerate() {
        // fetch the corresponding chrom index for this VCF file (they are not guaranteed to match)
        let vcf_header: &bcf::header::HeaderView = vcf_reader.header();
        let chrom_index: u32 = vcf_header.name2rid(region.get_chrom().as_bytes())?;

        // first make sure we find the sample in this file
        let mut lookup_index: Option<usize> = None;
        for (sample_index, &vcf_sample) in vcf_header.samples().iter().enumerate() {
            let vcf_sample_string: String = std::str::from_utf8(vcf_sample).unwrap().to_string();
            if vcf_sample_string == sample_name {
                lookup_index = Some(sample_index);
                break;
            }
        }
        match lookup_index {
            Some(index) => {
                sample_indices.push(index);
            },
            None => {
                bail!("Sample name {:?} was not found in VCF: {:?}", sample_name, vcf_paths[vcf_index]);
            }
        };
    
        // fetch our position in the VCF file
        match vcf_reader.fetch(chrom_index, region.get_start(), Some(region.get_end())) {
            Ok(()) => {
                // we have entries, so get the first one and queue it
                let mut vcf_iter = vcf_reader.records().peekable();
                let first_entry = vcf_iter.peek();
                if let Some(record_result) = first_entry {
                    let record: &rust_htslib::bcf::Record = match record_result {
                        Ok(r) => r,
                        // we have to convert to an owned error here, and the htslib errors are not cloneable
                        Err(e) => return Err(Box::new(SimpleError::from(e)))
                    };
                    let position: i64 = record.pos();
                    variant_queue.push(vcf_index, (Reverse(position), Reverse(vcf_index)));
                };

                // even if the iterator is empty, we push it so things are lined up correctly
                vcf_iterators.push(vcf_iter);
            },
            Err(_) => {
                // this usually happens when there are no entries for the chromosome
                vcf_iterators.push(vcf_reader.records().peekable());
            }
        };
    }

    // parse all the records and convert them into our format
    let mut variants: Vec<Variant> = Vec::<Variant>::with_capacity(region.get_num_variants());
    let mut hom_variants: Vec<Variant> = vec![];
    let mut previous_het_end: usize = 0;
    
    while !variant_queue.is_empty() {
        // get the source of the next variant to process and the sample index in that VCF file
        let (pop_index, pop_priority) = variant_queue.pop().unwrap();
        let sample_index: usize = sample_indices[pop_index];

        // process this variant
        let record_result = vcf_iterators[pop_index].next().unwrap();
        let record = record_result?;
        
        let position: i64 = record.pos();
        assert_eq!(position, pop_priority.0.0); // sanity check that the variant matches our position priority
        if position < region.get_start() as i64 {
            // this can happen when you have very very long indels that span one of our breaks
            // we have already written though, so don't write it again
        } else {
            let include_variant = is_phasable_variant(&record, sample_index, region.get_min_quality(), is_hom_allowed)?;
            if include_variant {
                let variant_type = get_variant_type(&record)?;
                
                // TODO: ideally, this would be consolidated with our Zygosity code block, but we need index_alleles further on
                //       possible solution is to make Zygosity types have u8 values tied to them
                //       low priority: not a major slowdown at this time
                // get the genotypes
                let all_genotypes = record.genotypes()?;
                let genotype = all_genotypes.get(sample_index);
                assert!(genotype.len() <= 2);

                // we don't really expect more than 255 alleles, but make sure we panic if that *does* happen
                let mut index_allele0: u8 = match genotype[0] {
                    GenotypeAllele::Unphased(at) => at.try_into().unwrap(),
                    GenotypeAllele::Phased(at) => at.try_into().unwrap(),
                    //TODO: ignore these for now, not sure how to handle it?
                    GenotypeAllele::UnphasedMissing => panic!("Should not happen"),
                    GenotypeAllele::PhasedMissing => panic!("Should not happen")
                };

                let mut index_allele1: u8 = if genotype.len() == 1 {
                    // TRGT can generate single-haplotype results, in this instance just copy index_allele0 and pretend it was homozygous
                    index_allele0
                } else {
                    match genotype[1] {
                        GenotypeAllele::Unphased(at) => at.try_into().unwrap(),
                        GenotypeAllele::Phased(at) => at.try_into().unwrap(),
                        //TODO: ignore these for now, not sure how to handle it?
                        GenotypeAllele::UnphasedMissing => panic!("Should not happen"),
                        GenotypeAllele::PhasedMissing => panic!("Should not happen")
                    }
                };

                // in merged VCF files, they are not always ordered, check for that here and swap if they are out of order
                // technically, this can happen in any VCF I supposed
                // are there any concerns to this swap? we already do many output swaps with the phasing so seems like no
                if index_allele0 > index_allele1 {
                    std::mem::swap(&mut index_allele0, &mut index_allele1);
                }

                // special case for homozygous
                let is_homozygous = index_allele0 == index_allele1;
                assert!(!is_homozygous || is_hom_allowed);
                if is_homozygous {
                    // this forces our homozygous variants to load as if they were heterozygous (e.g. with a reference allele)
                    // this is fine because we are not phasing them, just wanting to use the sequence
                    // TODO: can we remove this hackery? I could see it being an issue later
                    index_allele0 = 0;
                }
                
                let all_alleles = record.alleles();
                let ref_len: usize = all_alleles[0].len();
                let allele0: Vec<u8> = all_alleles[index_allele0 as usize].to_vec();
                let allele1: Vec<u8> = all_alleles[index_allele1 as usize].to_vec();

                let new_variant_result = match variant_type {
                    VariantType::Snv => {
                        Variant::new_snv(
                            pop_index, position, allele0, allele1, index_allele0, index_allele1
                        )
                    },
                    VariantType::Deletion => {
                        Variant::new_deletion(
                            pop_index, position, ref_len,
                            allele0, allele1, 
                            index_allele0, index_allele1
                        )
                    },
                    VariantType::Insertion => {
                        Variant::new_insertion(
                            pop_index, position, allele0, allele1, index_allele0, index_allele1
                        )
                    },
                    VariantType::Indel => {
                        Variant::new_indel(
                            pop_index, position, ref_len,
                            allele0, allele1, 
                            index_allele0, index_allele1
                        )
                    },
                    VariantType::SvDeletion => {
                        Variant::new_sv_deletion(
                            pop_index, position, ref_len,
                            allele0, allele1, 
                            index_allele0, index_allele1
                        )
                    },
                    VariantType::SvInsertion => {
                        Variant::new_sv_insertion(
                            pop_index, position, ref_len, allele0, allele1, index_allele0, index_allele1
                        )
                    },
                    VariantType::TandemRepeat => {
                        Variant::new_tandem_repeat(
                            pop_index, position, ref_len,
                            allele0, allele1, 
                            index_allele0, index_allele1
                        )
                    }
                    VariantType::SvDuplication |
                    VariantType::SvInversion |
                    VariantType::SvBreakend |
                    VariantType::Unknown => {
                        // panic here because we shouldn't allow these types unless we implement the variant
                        panic!("no impl for {variant_type:?}");
                    }
                };

                let mut new_variant = match new_variant_result {
                    Ok(nv) => nv,
                    Err(e) => {
                        bail!("Error processing variant in VCF#{pop_index} at {}:{} : {e}", region.get_chrom(), position+1);
                    }
                };

                if reference_buffer > 0 && !is_homozygous {
                    // we have a reference genome and a desire buffer, extend our alleles
                    let mut ref_prefix_start: usize = (position as usize).saturating_sub(reference_buffer);
                    let ref_postfix_start: usize = position as usize + ref_len;
                    
                    // we used to have an assertion here with the plan to remove it eventually, but it turns out
                    // that users do crazy things, so we need to convert it to a full Error
                    let ref_sequence = reference_genome.get_slice(region.get_chrom(), position as usize, ref_postfix_start);
                    if all_alleles[0] != ref_sequence {
                        // check there are IUPAC in the reference 
                        let iupac_filter_ref_sequence: Vec<u8> = ref_sequence.iter()
                            .map(|&c| {
                                match c {
                                    b'A' | b'C' | b'G' | b'T' => c,
                                    _ => b'N'
                                }
                            }).collect();
                        if all_alleles[0] != iupac_filter_ref_sequence {
                            bail!(
                                "Reference mismatch error: variant at {}:{} has REF allele = \"{}\", but reference genome has \"{}\".",
                                region.get_chrom(), position+1, 
                                // we don't want to panic in the middle of this, so use a safe unwrapper with default
                                std::str::from_utf8(all_alleles[0]).unwrap_or("utf8 decode error"), 
                                std::str::from_utf8(ref_sequence).unwrap_or("utf8 decode error")
                            );
                        } else {
                            // this is a case where a non-ACGT was replaced with N, so allow it
                            // debug!("IUPAC match");
                        }
                    }

                    // check if this variant is too close to the previous variant
                    if ref_prefix_start < previous_het_end {
                        // for the previous variant, we need to truncate it, possibly all the way down
                        if let Some(v) = variants.last_mut() {
                            let current_end: usize = v.position() as usize + v.get_ref_len() + v.get_postfix_len();
                            let truncate_length: usize = (current_end - position as usize).min(v.get_postfix_len());
                            v.truncate_reference_postfix(truncate_length);
                        } else {
                            panic!("This should not happen with our checks.");
                        }

                        // for this variant, set the new end to either the previous het end OR the position, whichever is lower
                        // previous het end can be lower if you have overlapping indels, not much we can do to fix that without global realignment
                        ref_prefix_start = previous_het_end.min(position as usize);
                    }

                    // add the prefix
                    let prefix: &[u8] = reference_genome.get_slice(region.get_chrom(), ref_prefix_start, position as usize);
                    new_variant.add_reference_prefix(prefix);

                    // add the postfix
                    let postfix: &[u8] = reference_genome.get_slice(region.get_chrom(), ref_postfix_start, ref_postfix_start + reference_buffer);
                    new_variant.add_reference_postfix(postfix);

                    // update the previous position to match the position + ref length
                    previous_het_end = position as usize + ref_len;
                }

                if is_homozygous {
                    hom_variants.push(new_variant);
                } else {
                    variants.push(new_variant);
                }
            }
        }

        // requeue from the one we popped from
        let next_entry = vcf_iterators[pop_index].peek();
        if let Some(record_result) = next_entry {
            let record: &rust_htslib::bcf::Record = match record_result {
                Ok(r) => r,
                // we have to convert to an owned error here, and the htslib errors are not cloneable
                Err(e) => return Err(Box::new(SimpleError::from(e)))
            };
            let position: i64 = record.pos();
            variant_queue.push(pop_index, (Reverse(position), Reverse(pop_index)));
        };
    }

    // sanity check that we found the same number of things
    assert_eq!(variants.len(), region.get_num_variants());
    Ok((variants, hom_variants))
}

/// A result for a phasing algorithm, assumes diploid solution currently.
pub struct PhaseResult {
    /// The phase block defining the problem space.
    pub phase_block: PhaseBlock,
    /// The variants contained in the phase block.
    pub variants: Vec<Variant>,
    /// The first haplotype in the solution.
    pub haplotype_1: Vec<AlleleType>,
    /// The second haplotype in the solution.
    pub haplotype_2: Vec<AlleleType>,
    /// Store the phase block ID of the variant, None if the variant is unphased
    pub block_ids: Vec<Option<usize>>,
    /// Stores all non-empty sub-blocks
    pub sub_phase_blocks: Vec<PhaseBlock>,
    /// Optional read statistics
    pub read_statistics: Option<ReadStats>,
    /// Optional statistics from the problem
    pub statistics: Option<PhaseStats>
}

/// Calculates the span count for each juncture in the solution, ignoring homozygous variants and unassigned alleles.
/// # Arguments
/// * `read_segments` - the read segments used to solve the problem
/// * `haplotype_1` - the first haplotype in the solution
/// * `haplotype_2` - the second haplotype in the solution
fn get_solution_span_counts(
    read_segments: &IntervalTree<usize, ReadSegment>,
    haplotype_1: &[AlleleType],
    haplotype_2: &[AlleleType]
) -> Vec<usize> {
    // this will store the total spanning reads ignoring any homozygous or unsolved variants
    assert_eq!(haplotype_1.len(), haplotype_2.len());

    // there is one less connection than the total number of variants
    let mut total_span_counts: Vec<usize> = vec![0; haplotype_1.len() - 1];
    
    // iterate over all the read segments
    for rs_interval in read_segments.find(0..usize::MAX) {
        let rs = rs_interval.data();

        // the range returns [first_allele, last_allele+1), we need to basically remove the +1 here since we're talking junctures
        let mut juncture_range = rs.region().clone();
        juncture_range.end -= 1;

        // if any of the head variants were converted to homozygous, do not include because they don't provide spanning evidence anymore
        while juncture_range.start < juncture_range.end && 
            haplotype_1[juncture_range.start] == haplotype_2[juncture_range.start] {
                juncture_range.start += 1;
        }

        // if any of the tail variants were converted to homozygous, do not include because they don't provide spanning evidence anymore
        while juncture_range.start < juncture_range.end &&
            haplotype_1[juncture_range.end] == haplotype_2[juncture_range.end] {
                juncture_range.end -= 1;
        }
        
        // the range has been truncated upstream so we're only looking at junctures
        for tpc in total_span_counts[juncture_range].iter_mut() {
            *tpc += 1;
        }
    }
    
    total_span_counts
}

/// Core structure of phasing that can be run on a single processor to solve a phase block.
/// This method is designed to perform data loading and then run an algorithm from another module to solve the block.
/// See `astar_phaser::astar_solver(...)` for an example solver implementation.
/// # Arguments
/// * `phase_problem` - the problem definition, primarily defines coordinates of the phase block we want to solve
/// * `vcf_paths` - the VCF files to load variants from, must be zipped and indexed
/// * `sample_name` - the sample name inside the VCF files
/// * `bam_paths` - the BAM files to load read observations from, must be indexed
/// * `reference_genome` - optional, the reference genome 
/// * `reference_buffer` - the number of nearby bases to try to use for local realignment
/// * `min_matched_alleles` - the minimum number of matched alleles required to include a read
/// * `min_mapq` - the minimum MAPQ to include a read
/// * `astar_config` - the configuration for the A* algorithm
/// * `global_config` - the global realignment configuration; if None, then only local realignment is used
pub fn solve_block(
    phase_problem: &PhaseBlock, vcf_paths: &[PathBuf], bam_paths: &[PathBuf], 
    reference_genome: &ReferenceGenome, phasing_config: &PhasingConfig
) -> Result<(PhaseResult, HaplotagResult), Box<dyn std::error::Error>> {
    debug!("Solving problem: {:?}", phase_problem);
    // unwrap the phasing configuration
    let reference_buffer: usize = phasing_config.reference_buffer;
    let min_matched_alleles: usize = phasing_config.min_matched_alleles;
    let min_mapq: u8 = phasing_config.min_mapq;
    let astar_config: &AstarConfig = &phasing_config.astar_config;
    let global_config: Option<&GlobalRealignmentConfig> = phasing_config.global_realignment_config.as_ref();

    // short circuit for "empty" problems
    if phase_problem.get_num_variants() == 0 {
        // this should only happen for chromosomes with no het alleles
        assert!(phase_problem.get_start() == 0);
        assert!(phase_problem.get_end() == 0);
        let empty_result = PhaseResult {
            phase_block: phase_problem.clone(),
            variants: vec![],
            haplotype_1: vec![],
            haplotype_2: vec![],
            block_ids: vec![],
            sub_phase_blocks: vec![],
            read_statistics: None,
            statistics: None
        };
        let empty_haplotag_result: HaplotagResult = HaplotagResult { 
            phase_block: phase_problem.clone(),
            reads: Default::default() 
        };
        return Ok((empty_result, empty_haplotag_result));
    }

    // homs are only used if global realignment is being attempted
    let load_homs: bool = global_config.is_some();

    // lets extract the variants we care about from the vcf
    let (mut variant_calls, mut hom_calls): (Vec<Variant>, Vec<Variant>) = load_variant_calls(
        phase_problem, 
        vcf_paths, 
        reference_genome, reference_buffer, 
        load_homs
    )?;
    assert_eq!(variant_calls.len(), phase_problem.get_num_variants());

    // go through all the loaded variants, including homs, and pull out the TandemRepeat coordinates that have been loaded
    let mut tr_segments: IntervalTree<i64, usize> = Default::default();
    for variant in variant_calls.iter().chain(hom_calls.iter()) {
        if variant.get_type() == VariantType::TandemRepeat {
            let start: i64 = variant.position();
            let ref_len: usize = variant.get_ref_len();
            let end: i64 = start + ref_len as i64;
            // prior to TRGT v1.0.0: sometimes TRGT inserts the base before which will match an insertion and sometimes it wouldn't
            // now there is always an anchor base, so we can safely just use start..end from the VCF directly
            tr_segments.insert(start..end, 0);
        }
    }

    // now mark all variants contained by the STRs as ignored
    for variant in variant_calls.iter_mut() {
        if variant.get_type() != VariantType::TandemRepeat {
            let start: i64 = variant.position();
            let ref_len: usize = variant.get_ref_len();
            let end: i64 = start + ref_len as i64;
            let var_interval = start..end;
            let mut is_contained: bool = false;
            for segment in tr_segments.find(var_interval) {
                let seg_start = segment.interval().start;
                let seg_end = segment.interval().end;
                if seg_start <= start && seg_end >= end {
                    // this segment fully contains the variant
                    is_contained = true;
                    break;
                }
            }

            if is_contained {
                // we need to mark this one as unphaseable
                variant.set_ignored(IgnoredVariantReason::TandemRepeatOverlap);
                debug!("Set ignored het: {:?}", variant);
            }
        }
    }

    // mark any homs we are going to ignore as well
    for variant in hom_calls.iter_mut() {
        if variant.get_type() != VariantType::TandemRepeat {
            let start: i64 = variant.position();
            let ref_len: usize = variant.get_ref_len();
            let end: i64 = start + ref_len as i64;
            let var_interval = start..end;
            let mut is_contained: bool = false;
            for segment in tr_segments.find(var_interval) {
                let seg_start = segment.interval().start;
                let seg_end = segment.interval().end;
                if seg_start <= start && seg_end >= end {
                    // this segment fully contains the variant
                    is_contained = true;
                    break;
                }
            }

            if is_contained {
                // we need to mark this one as unphaseable
                variant.set_ignored(IgnoredVariantReason::TandemRepeatOverlap);
                debug!("Set ignored hom: {:?}", variant);
            }
        }
    }

    // reads that meet our criteria for use in phasing
    let read_segments: IntervalTree<usize, ReadSegment>;
    // reads that meet our criteria EXCEPT for the minimum number of alleles is just > 0, so they can potentially be phased after we solve
    let phasable_segments: IntervalTree<usize, ReadSegment>;
    let read_stats: ReadStats;

    if let Some(grc) = global_config {
        
        // we are attempting global re-alignments
        (read_segments, phasable_segments, read_stats) = read_parsing::load_full_read_segments(
            phase_problem, bam_paths, &variant_calls, &hom_calls,
            reference_genome, min_matched_alleles, min_mapq,
            grc
        )?;
    } else {
        // we are doing local re-alignments only
        (read_segments, phasable_segments, read_stats) = read_parsing::load_read_segments(
            phase_problem, bam_paths, reference_genome.filename(),
            &variant_calls, min_matched_alleles, min_mapq
        )?;
    }

    // TODO: here we should check if there are any variants where no alleles are set; if so, we should mark them as ignored
    //   furthermore, this may require clearing those alleles, but we can handle that later

    // count the assigned variants for each read segment
    let assigned_counts: Vec<usize> = count_assigned_variants(variant_calls.len(), &read_segments);
    for (i, &count) in assigned_counts.iter().enumerate() {
        if count == 0 {
            // this variant is not assigned to any reads
            variant_calls[i].set_ignored(IgnoredVariantReason::NoReadsAssigned);
            debug!("Set ignored variant, no reads assigned: {:?}", variant_calls[i]);
        }
    }

    // read segment debugging
    for (i, seg) in read_segments.find(0..usize::MAX).enumerate() {
        trace!("read segment #{} => {:?}", i, seg);
    }

    // okay final phase is to solve some algorithm given those read segments
    let reorder_variants = phasing_config.optimize_variant_order;
    let astar_result: astar_phaser::AstarResult = if reorder_variants {
        // optimize the variant order
        let reordered_problem = optimize_variant_order(&variant_calls[..], &read_segments)?;

        // run our solver with the reordered problem
        let reordered_astar_result = astar_phaser::astar_solver(
            phase_problem, // the phase block summary information, this should be okay here
            reordered_problem.variants(),
            reordered_problem.read_segments(),
            astar_config
        );

        // translate the AstarResult to match the original variant order
        reordered_problem.translate_solution(&reordered_astar_result)
    } else {
        astar_phaser::astar_solver(
            phase_problem, &variant_calls[..], &read_segments, astar_config
        )
    };

    debug!("Building phase block IDs...");
    let block_tags = if phasing_config.use_linear_block_algorithm {
        // old linear block algorithm
        build_linear_phase_tags(
            &variant_calls, &read_segments, &astar_result
        )
    } else {
        // non-linear block algorithm, with a configurable minimum number of connecting reads
        let min_connecting_reads = phasing_config.min_connecting_reads;
        get_phase_block_ids(
            &variant_calls, &read_segments, &astar_result, min_connecting_reads
        )
    };
    debug!("block_ids: {:?}", block_tags);

    // we have the variants grouped into tags (not necessarily final tag value), now build the phase blocks
    // this function does not require linear blocks, and it works for both methods
    let sub_phase_blocks = build_sub_phase_blocks(
        phase_problem, &variant_calls, &block_tags
    );

    // last step is to haplotag the reads we loaded
    let mut haplotagged_reads: HashMap<String, (usize, usize)> = haplotag_reads(
        read_segments, &astar_result.haplotype_1, &astar_result.haplotype_2,
        &block_tags
    );

    // also haplotype the extra reads and add them to our result
    let phasable_haplotagged_reads: HashMap<String, (usize, usize)> = haplotag_reads(
        phasable_segments, &astar_result.haplotype_1, &astar_result.haplotype_2,
        &block_tags
    );

    // we could just extend here, but we want to sanity check our keys don't overlap at all
    for (k, v) in phasable_haplotagged_reads.into_iter() {
        // we should never have the same read name in both hashmaps
        assert!(!haplotagged_reads.contains_key(&k));
        haplotagged_reads.insert(k, v);
    }

    let haplotag_result: HaplotagResult = HaplotagResult {
        phase_block: phase_problem.clone(),
        reads: haplotagged_reads
    };

    // save all our results here
    let phase_result: PhaseResult = PhaseResult {
        phase_block: phase_problem.clone(),
        variants: variant_calls,
        haplotype_1: astar_result.haplotype_1,
        haplotype_2: astar_result.haplotype_2,
        block_ids: block_tags,
        sub_phase_blocks,
        read_statistics: Some(read_stats),
        statistics: Some(astar_result.statistics)
    };
    Ok((phase_result, haplotag_result))
}

/// Returns a vector of the number of assigned variants for each variant in the block.
/// # Arguments
/// * `num_variants` - the number of variants in the block
/// * `read_segments` - the read segments to count the assigned variants for
/// # Returns
/// A vector of the number of assigned variants for each variant in the block
fn count_assigned_variants(num_variants: usize, read_segments: &IntervalTree<usize, ReadSegment>) -> Vec<usize> {
    let mut assigned_counts: Vec<usize> = vec![0; num_variants];
    for rs_interval in read_segments.find(0..usize::MAX) {
        // get the segment and its region
        let rs = rs_interval.data();
        let region = rs.region();
        #[allow(clippy::needless_range_loop)] // much easier to read this way, ignore clippy
        for var_index in region.start..region.end {
            // variants are only assigned if they are definitively reference or alternate
            let allele = rs.allele(var_index);
            if allele == AlleleType::Reference || allele == AlleleType::Alternate {
                assigned_counts[var_index] += 1;
            }
        }
    }
    assigned_counts
}

/// This function generates a dummy result for a block we are not phasing.
/// It is boilerplate for an unsolved block, either because there is only one variant we do not wish to phase OR because there are no reads but there are variants.
/// In either case, we will mark all the variant alleles as 2/2 to indicate unphased.
/// # Arguments
/// * `phase_problem` - the problem definition, primarily defines coordinates of the phase block we want to solve
pub fn create_unphased_result(phase_problem: &PhaseBlock) -> (PhaseResult, HaplotagResult) {
    debug!("Generating unphased result for block: {phase_problem:?}");
    
    // generate a vec with the appropriate number of variants from each VCF
    let num_variants = phase_problem.get_num_variants();
    let mut variant_calls: Vec<Variant> = vec![];
    for (vcf_index, &vi_count) in phase_problem.vcf_index_counts().iter().enumerate() {
        let dummy_variant = Variant::new_snv(
            vcf_index,
            phase_problem.get_start() as i64,
            vec![0],
            vec![1],
            0,
            1
        ).unwrap();
        variant_calls.extend(std::iter::repeat_n(dummy_variant, vi_count));
    }

    assert_eq!(variant_calls.len(), num_variants);
    
    // now we can make our dummy results - 0/0 get left as an unphased het
    // originally, this was 2/2, but that's a sentinel for TR_OVERLAP; we should make that more obvious in the future
    let phase_result: PhaseResult = PhaseResult {
        phase_block: phase_problem.clone(),
        variants: variant_calls,
        haplotype_1: vec![AlleleType::Reference; num_variants],
        haplotype_2: vec![AlleleType::Reference; num_variants],
        block_ids: vec![None; num_variants],
        sub_phase_blocks: vec![], // empty because this is not getting treated as a block
        read_statistics: None,
        statistics: None
    };
    let haplotag_result: HaplotagResult = HaplotagResult { 
        phase_block: phase_problem.clone(), 
        reads: Default::default() // no haplotagging in this mode, so give back an empty map
    };
    (phase_result, haplotag_result)
}

/// Builds the linear phase tags for the variants in the block.
/// This is the original linear approach where we assume blocks are contiguous and non-overlapping.
/// # Arguments
/// * `variant_calls` - the variant calls to build the phase tags for
/// * `read_segments` - the read segments to use to build the phase tags
/// * `astar_result` - the AstarResult to use to build the phase tags
/// # Returns
/// A vector of the phase tags for the variants in the block
fn build_linear_phase_tags(
    variant_calls: &[Variant],
    read_segments: &IntervalTree<usize, ReadSegment>,
    astar_result: &AstarResult
) -> Vec<Option<usize>> {
    // get the total spanning counts after accounting for homozygous variants / missing alleles
    let total_span_counts: Vec<usize> = get_solution_span_counts(read_segments, &astar_result.haplotype_1, &astar_result.haplotype_2);
    debug!("total_span_counts: {:?}", total_span_counts);

    // if we end up having no reads spanning a juncture, it's time to split the block
    let block_split: Vec<bool> = total_span_counts.iter()
        .map(|&tc| tc == 0)
        .collect::<Vec<bool>>();

    debug!("Block split: {:?}", block_split);

    // now, figure out what the haplotag is for each variant
    let mut block_tags: Vec<Option<usize>> = vec![None; variant_calls.len()];
    let mut current_tag: usize = variant_calls[0].position() as usize;
    let haplotype_1 = &astar_result.haplotype_1;
    let haplotype_2 = &astar_result.haplotype_2;
    for (i, variant) in variant_calls.iter().enumerate() {
        // we will not tag any ignored variants or variants that are flagged as homozygous
        if !variant.is_ignored() && haplotype_1[i] != haplotype_2[i] {
            if i > 0 && block_split[i-1] {
                // this is a new block
                current_tag = variant.position() as usize;
            }

            // set the block tag
            block_tags[i] = Some(current_tag);
        }
    }
    debug!("Block tags: {:?}", block_tags);
    block_tags
}

/// Builds the phase block IDs for the variants in the block.
/// This is a new approach where we can have blocks that overlap due to things like splicing.
/// # Arguments
/// * `variants` - the variant calls to build the phase block IDs for
/// * `read_segments` - the read segments to use to build the phase block IDs
/// * `astar_result` - the AstarResult to use to build the phase block IDs
/// * `min_connecting_reads` - the minimum number of reads that must connect a variant to a block to be included in the block
/// # Returns
/// A vector of the phase block IDs for the variants in the block
/// # Panics
/// * If `min_connecting_reads` is less than 1
fn get_phase_block_ids(
    variants: &[Variant],
    read_segments: &IntervalTree<usize, ReadSegment>,
    astar_result: &AstarResult,
    min_connecting_reads: u64
) -> Vec<Option<usize>> {
    assert!(min_connecting_reads > 0);

    // get the haplotypes from the AstarResult
    let haplotype_1 = &astar_result.haplotype_1;
    let haplotype_2 = &astar_result.haplotype_2;

    // sanity checks
    assert_eq!(variants.len(), haplotype_1.len());
    assert_eq!(variants.len(), haplotype_2.len());

    // initialize the block ids to none
    let num_variants: usize = variants.len();
    let mut block_ids: Vec<Option<usize>> = vec![None; num_variants];

    for i in 0..num_variants {
        if haplotype_1[i] == haplotype_2[i] || 
            block_ids[i].is_some() {
            // this is a homozygous variant, so we don't need to check any further OR
            // we already have a block id for this variant, so we can skip to the next variant
            continue;
        }

        // initialize our block ID to the position of this variant
        let block_id = variants[i].position() as usize;
        let mut variants_in_block: HashSet<usize> = HashSet::default();
        let mut variants_to_explore = BTreeSet::default();
        variants_to_explore.insert(i);

        while let Some(variant_index) = variants_to_explore.pop_first() {
            if variants_in_block.contains(&variant_index) {
                // this variant is already in our block, which means we already analyzed it
                // do a quick sanity check and move on
                assert_eq!(block_ids[variant_index], Some(block_id));
                continue;
            }

            // make sure this variant is not ignored
            assert!(!variants[variant_index].is_ignored());
            assert_ne!(haplotype_1[variant_index], haplotype_2[variant_index]);

            // add this variant to the block
            variants_in_block.insert(variant_index);
            block_ids[variant_index] = Some(block_id);

            debug!("Adding variant {} to block {}", variant_index, block_id);
            let mut connecting_counts: BTreeMap<usize, u64> = BTreeMap::default();

            // now go through all the read segments that contain this variant
            for rs_interval in read_segments.find(variant_index..variant_index+1) {
                // get the read segment
                let rs = rs_interval.data();
                let rs_range = rs.region();

                // check that the variant we popped is a valid allele for this read segment
                let i_allele = rs.allele(variant_index);
                if i_allele == AlleleType::Ambiguous || i_allele == AlleleType::NoOverlap {
                    // this read segment does not contain a valid allele for this variant, so we can skip it
                    continue;
                }

                // iterate over the other valid alleles in the read segment
                for j in rs_range.start..rs_range.end {
                    // get the allele at this index
                    let allele = rs.allele(j);
                    if allele == AlleleType::Ambiguous || allele == AlleleType::NoOverlap || // read is ambiguous or no overlap at this allele
                        haplotype_1[j] == haplotype_2[j] || // this allele is homozygous in our solution
                        variants_in_block.contains(&j) { // this allele is already in our block, do not double count
                        continue;
                    } else {
                        // add the variant to the list of variants to explore
                        // variants_to_explore.insert(j);
                        *connecting_counts.entry(j).or_insert(0) += 1;
                    }
                }
            }

            // now, check if any of the variants have enough connecting reads to be included in the block
            for (variant_index, count) in connecting_counts.iter() {
                if *count >= min_connecting_reads {
                    variants_to_explore.insert(*variant_index);
                }
            }
        }
    }

    // return the block ids
    block_ids
}

fn build_sub_phase_blocks(
    phase_problem: &PhaseBlock,
    variant_calls: &[Variant],
    block_tags: &[Option<usize>]
) -> Vec<PhaseBlock> {
    // initialize our sub-phase blocks, we want these to be sorted by the block index so use a BTreeMap
    let mut sub_phase_blocks: BTreeMap<usize, PhaseBlock> = BTreeMap::default();

    // now, iterate through the variant calls and build the sub-phase blocks
    for (variant, block_tag) in variant_calls.iter().zip(block_tags.iter()) {
        if let Some(block_index) = block_tag {
            let block = sub_phase_blocks.entry(*block_index)
                .or_insert(PhaseBlock::new(
                    phase_problem.get_block_index(),
                    phase_problem.get_chrom().to_string(),
                    phase_problem.get_chrom_index(),
                    phase_problem.get_min_quality(),
                    phase_problem.sample_name().to_string(),
                    phase_problem.vcf_index_counts().len()
                ));

            // add the variant to the current block
            block.add_locus_variant(
                phase_problem.get_chrom(), variant.position() as u64, variant.get_vcf_index()
            );
        }
    }

    // convert the BTreeMap to a Vec of the values, we can drop the keys now
    sub_phase_blocks.into_values().collect()
}

/// Stores all information for a haplotag result
#[derive(Debug)]
pub struct HaplotagResult {
    /// The phase block defining the problem space.
    pub phase_block: PhaseBlock,
    /// Indexes reads by name and returns a tuple (phase block ID, haplotag)
    pub reads: HashMap<String, (usize, usize)>
}

/// Returns the tagging results for the reads in a HashMap.
/// Values are tuples (phase block ID (0-based), haplotag value (0 or 1))
/// # Arguments
/// * `read_segments` - the reads to tag
/// * `haplotype_1` - the first haplotype
/// * `haplotype_2` - the second haplotype
/// * `block_tags` - tags for the blocks based on the variant IDs
/// # Panics
/// * if the haplotypes are not the same length as each other, a read segment, and/or the variant calls
/// * the block breaks is not 1 less length than the variant calls
pub fn haplotag_reads(
    read_segments: IntervalTree<usize, ReadSegment>, 
    haplotype_1: &[AlleleType], haplotype_2: &[AlleleType], block_tags: &[Option<usize>]
) -> HashMap<String, (usize, usize)> {
    // now do the tagging
    let mut haplotagged_reads: HashMap<String, (usize, usize)> = Default::default();
    for rs_interval in read_segments.find(0..usize::MAX) {
        // first, see if we can resolve it to a haplotype
        let rs: &ReadSegment = rs_interval.data();
        let a1_score: u64 = rs.score_haplotype(haplotype_1);
        let a2_score: u64 = rs.score_haplotype(haplotype_2);
        let haplotag: usize = match a1_score.cmp(&a2_score) {
            Ordering::Less => 0,
            Ordering::Greater => 1,
            Ordering::Equal => 2
        };

        if haplotag != 2 {
            // we can resolve to a haplotype, now get the phase block index
            // find the first resolved variant in our read segment
            let region = rs.region();
            let mut first_variant: usize = region.start;

            // TODO: it's possible that a read hits multiple blocks; this approach tags by the first variant with a tag
            //       this may not be the optimal method, but we'll leave it for now

            // the variant is untagged OR the variant is not resolved (which can happen sometimes)
            while block_tags[first_variant].is_none() || rs.allele(first_variant) >= AlleleType::Ambiguous {
                first_variant += 1;
            }
            assert!(first_variant < region.end); // this should never happen, but we'll assert it just in case

            // assuming the phase block is set, we can tag the read now
            if let Some(phase_block) = block_tags[first_variant] {
                // finally, just get the read name and make sure we haven't somehow already marked this one
                let read_name: String = rs.read_name().to_string();
                assert!(!haplotagged_reads.contains_key(&read_name));
                haplotagged_reads.insert(read_name, (phase_block, haplotag));
            }
        }
    }

    haplotagged_reads
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_solution_span_counts() {
        let haplotype_1: Vec<AlleleType> = vec![0, 1, 1, 0, 0, 0].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect();
        let haplotype_2: Vec<AlleleType> = vec![1, 1, 1, 1, 0, 1].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect();
        let test_reads = vec![
            ReadSegment::new("r1".to_string(), vec![0, 0, 0, 0, 0, 0].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect(), vec![1, 1, 1, 1, 1, 1]), // adds 1 to everything
            ReadSegment::new("r2".to_string(), vec![3, 3, 3, 1, 1, 3].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect(), vec![0, 0, 0, 1, 1, 0]), // one allele is a hom, so this does nothing
            ReadSegment::new("r3".to_string(), vec![1, 1, 1, 1, 3, 3].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect(), vec![1, 1, 1, 1, 0, 0]), // adds 1 to first 3
            ReadSegment::new("r4".to_string(), vec![3, 1, 1, 1, 1, 1].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect(), vec![0, 1, 1, 1, 1, 1]), // adds 1 to last 2
        ];
        let mut read_segments: IntervalTree<usize, ReadSegment> = Default::default();
        for rs in test_reads.into_iter() {
            let rs_range = rs.region().clone();
            read_segments.insert(rs_range, rs);
        }
        let expected_result: Vec<usize> = vec![2, 2, 2, 2, 2];

        let result = get_solution_span_counts(&read_segments, &haplotype_1, &haplotype_2);
        assert_eq!(expected_result, result);
    }

    #[test]
    fn test_haplotag_reads() {
        let haplotype_1: Vec<AlleleType> = vec![0, 0, 0, 0, 0, 0].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect();
        let haplotype_2: Vec<AlleleType> = vec![1, 1, 1, 1, 1, 1].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect();
        let block_tags = vec![Some(0), Some(0), Some(0), Some(3), Some(3), Some(5)];
        let test_reads = vec![
            ReadSegment::new("r1".to_string(), vec![0, 0, 0, 0, 0, 0].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect(), vec![1, 1, 1, 1, 1, 1]),
            ReadSegment::new("r2".to_string(), vec![2, 2, 2, 1, 1, 2].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect(), vec![0, 0, 0, 1, 1, 0]),
            ReadSegment::new("r3".to_string(), vec![2, 2, 2, 1, 0, 2].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect(), vec![0, 0, 0, 1, 1, 0]),
            ReadSegment::new("r4".to_string(), vec![2, 2, 2, 1, 0, 1].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect(), vec![0, 0, 0, 1, 1, 1]),
            ReadSegment::new("r5".to_string(), vec![2, 2, 2, 1, 0, 2].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect(), vec![0, 0, 0, 2, 1, 0]),
        ];

        let mut read_segments: IntervalTree<usize, ReadSegment> = Default::default();
        for rs in test_reads.into_iter() {
            let rs_range = rs.region().clone();
            read_segments.insert(rs_range, rs);
        }

        let haplotag_result = haplotag_reads(read_segments, &haplotype_1, &haplotype_2, &block_tags);

        // simple-ish cases
        assert_eq!(haplotag_result.get("r1").unwrap(), &(0, 0)); // exact match to haplotype 1
        assert_eq!(haplotag_result.get("r2").unwrap(), &(3, 1)); // exact, but incomplete, match to haplotype 2
        assert!(!haplotag_result.contains_key("r3")); // equal to both, so unassigned
        assert_eq!(haplotag_result.get("r4").unwrap(), &(3, 1)); // spans two blocks and is inexact, but closer to hap 1 starting with block 3
        assert_eq!(haplotag_result.get("r2").unwrap(), &(3, 1)); // equal by alleles, but qual on allele 1 is higher
    }

    // tests that having a IUPAC code modified to "N" is okay
    #[test]
    fn test_iupac_conversion() {
        // create a single variant block for simplicity
        let mut phase_block = PhaseBlock::new(0, "chr1".to_string(), 0, 20, "sample_name".to_string(), 1);
        phase_block.add_locus_variant("chr1", 4, 0);

        // load it up
        let vcf_path = PathBuf::from("./test_data/iupac_test/small_variants.vcf.gz");
        let ref_fn = PathBuf::from("./test_data/iupac_test/test_reference.fa");
        let reference_genome = ReferenceGenome::from_fasta(&ref_fn).unwrap();
        let reference_buffer = 2;
        let is_hom_allowed = false;
        let (het_variants, hom_variants) = load_variant_calls(
            &phase_block,
            &[vcf_path],
            &reference_genome,
            reference_buffer,
            is_hom_allowed
        ).unwrap();

        // we should just have the one variant in the block with the translated IUPAC code
        let mut expected_variant = Variant::new_indel(0, 4, 4, b"ANGT".to_vec(), b"AT".to_vec(), 0, 1).unwrap();
        expected_variant.add_reference_prefix(b"GT");
        expected_variant.add_reference_postfix(b"AC");
        
        assert_eq!(het_variants, vec![
            expected_variant
        ]);
        assert!(hom_variants.is_empty());
    }

    /// Returns a generic AstarResult with the given length for the haplotypes.
    /// # Arguments
    /// * `length` - the length of the haplotypes
    /// # Returns
    /// A generic AstarResult with the given length for the haplotypes
    fn get_generic_astar_result(length: usize) -> AstarResult {
        AstarResult {
            haplotype_1: vec![AlleleType::Reference; length],
            haplotype_2: vec![AlleleType::Alternate; length],
            statistics: PhaseStats::astar_new(0, 0, 0, length as u64, length as u64, 0, 0),
        }
    }

    #[test]
    fn test_complete_phase_tags() {
        // test case 1: single continuous block - all variants connected
        let astar_result = get_generic_astar_result(3);
        let variants1 = vec![
            Variant::new_snv(0, 100, b"A".to_vec(), b"T".to_vec(), 0, 1).unwrap(),
            Variant::new_snv(0, 200, b"G".to_vec(), b"C".to_vec(), 0, 1).unwrap(),
            Variant::new_snv(0, 300, b"T".to_vec(), b"A".to_vec(), 0, 1).unwrap(),
        ];
        let test_reads = vec![
            // Read spanning all three variants
            ReadSegment::new("r1".to_string(),
                vec![0, 0, 0].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect(),
                vec![10, 10, 10]),
            ReadSegment::new("r2".to_string(),
                vec![1, 1, 1].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect(),
                vec![10, 10, 10]),
        ];
        let mut read_segments: IntervalTree<usize, ReadSegment> = Default::default();
        for rs in test_reads.into_iter() {
            let rs_range = rs.region().clone();
            read_segments.insert(rs_range, rs);
        }

        // All variants should be in the same block (tagged with first variant's position)
        let linear_result = build_linear_phase_tags(&variants1, &read_segments, &astar_result);
        assert_eq!(linear_result, vec![Some(100), Some(100), Some(100)]);
        let clump_result = get_phase_block_ids(&variants1, &read_segments, &astar_result, 1);
        assert_eq!(clump_result, vec![Some(100), Some(100), Some(100)]);

        // check that we can require more than 1 read to connect a variant to a block
        let clump_result = get_phase_block_ids(&variants1, &read_segments, &astar_result, 3);
        assert_eq!(clump_result, vec![Some(100), Some(200), Some(300)]);
    }

    #[test]
    fn test_interwoven_phase_tags() {
        // test case: two overlapping blocks that are separate
        let astar_result = get_generic_astar_result(3);
        let variants1 = vec![
            Variant::new_snv(0, 100, b"A".to_vec(), b"T".to_vec(), 0, 1).unwrap(),
            Variant::new_snv(0, 200, b"G".to_vec(), b"C".to_vec(), 0, 1).unwrap(),
            Variant::new_snv(0, 300, b"T".to_vec(), b"A".to_vec(), 0, 1).unwrap(),
        ];
        let test_reads = vec![
            // Read spanning all three variants
            ReadSegment::new("r1".to_string(),
                vec![0, 3, 0].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect(),
                vec![10, 0, 10]),
            ReadSegment::new("r2".to_string(),
                vec![3, 1, 3].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect(),
                vec![0, 10, 0]),
        ];
        let mut read_segments: IntervalTree<usize, ReadSegment> = Default::default();
        for rs in test_reads.into_iter() {
            let rs_range = rs.region().clone();
            read_segments.insert(rs_range, rs);
        }

        // linear result thinks they all go together, which is why it's not as good
        let linear_result = build_linear_phase_tags(&variants1, &read_segments, &astar_result);
        assert_eq!(linear_result, vec![Some(100), Some(100), Some(100)]);

        // this one finds the two clumps
        let clump_result = get_phase_block_ids(&variants1, &read_segments, &astar_result, 1);
        assert_eq!(clump_result, vec![Some(100), Some(200), Some(100)]);
    }

    #[test]
    fn test_homozygous_phase_tags() {
        // test case: first variant is homozygous, so it's not in the block and not tagged
        let astar_result = AstarResult {
            haplotype_1: vec![AlleleType::Reference, AlleleType::Reference, AlleleType::Reference],
            haplotype_2: vec![AlleleType::Reference, AlleleType::Alternate, AlleleType::Alternate],
            statistics: PhaseStats::astar_new(0, 0, 0, 3, 3, 0, 0),
        };
        let variants1 = vec![
            Variant::new_snv(0, 100, b"A".to_vec(), b"T".to_vec(), 0, 1).unwrap(),
            Variant::new_snv(0, 200, b"G".to_vec(), b"C".to_vec(), 0, 1).unwrap(),
            Variant::new_snv(0, 300, b"T".to_vec(), b"A".to_vec(), 0, 1).unwrap(),
        ];
        let test_reads = vec![
            // Read spanning all three variants
            ReadSegment::new("r1".to_string(),
                vec![0, 0, 0].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect(),
                vec![10, 10, 10]),
            ReadSegment::new("r2".to_string(),
                vec![0, 1, 1].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect(),
                vec![0, 10, 10]),
        ];
        let mut read_segments: IntervalTree<usize, ReadSegment> = Default::default();
        for rs in test_reads.into_iter() {
            let rs_range = rs.region().clone();
            read_segments.insert(rs_range, rs);
        }

        // linear result thinks they all go together, which is why it's not as good
        let linear_result = build_linear_phase_tags(&variants1, &read_segments, &astar_result);
        assert_eq!(linear_result, vec![None, Some(200), Some(200)]);

        // this one finds the two clumps
        let clump_result = get_phase_block_ids(&variants1, &read_segments, &astar_result, 1);
        assert_eq!(clump_result, vec![None, Some(200), Some(200)]);
    }
}