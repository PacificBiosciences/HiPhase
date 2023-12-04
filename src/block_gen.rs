
use crate::data_types::variants::{VariantType, Zygosity};

use log::{debug, trace, warn};
use priority_queue::PriorityQueue;
use rust_htslib::{bam,bcf,htslib};
use rust_htslib::bcf::record::GenotypeAllele;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use simple_error::{SimpleError, bail};
use std::cell::RefCell;
use std::cmp::Reverse;
use std::path::{Path, PathBuf};

/// Uses for phase block priority during multi iteration
type PhaseBlockPriority = Reverse<(u32, u64, u64)>;

/// Gets a list of sample names from a given VCF file
/// # Arguments
/// * `filename` - the VCF file to load
/// # Errors
/// * if the filename fails to load as a VCF
/// * if the sample name fails to parse from utf8
pub fn get_vcf_samples(filename: &Path) -> Result<Vec<String>, Box<dyn std::error::Error>> {
    use rust_htslib::bcf::Read;
    let vcf_reader: bcf::IndexedReader = bcf::IndexedReader::from_path(filename)?;
    let vcf_header: bcf::header::HeaderView = vcf_reader.header().clone();
    let mut sample_names = vec![];
    for sv in vcf_header.samples().iter() {
        let vcf_sample_string: String = std::str::from_utf8(sv)?.to_string();
        sample_names.push(vcf_sample_string);
    }
    Ok(sample_names)
}

/// Iterates through a collection of BAM files and finds the ones matching the given sample name
/// # Arguments
/// * `all_bam_files` - the full list of BAM files to parse
/// * `sample_name` - the name we want to match
/// * `reference_filename` - the path to the reference file
/// # Errors
/// * if BAM has no RG tags
/// * if an RG tag has no SM tag
/// * if multiple SM tags are detected
pub fn get_sample_bams(all_bam_files: &[PathBuf], sample_name: &str, reference_filename: &Path) -> Result<(Vec<PathBuf>, Vec<usize>), Box<dyn std::error::Error>> {
    use rust_htslib::bam::Read;
    let mut sample_bams: Vec<PathBuf> = vec![];
    let mut bam_indices: Vec<usize> = vec![];
    for (bam_index, bam_fn) in all_bam_files.iter().enumerate() {
        let bam_reader = {
            let mut b = bam::IndexedReader::from_path(bam_fn)?;
            b.set_reference(reference_filename)?;
            b
        };
        let bam_header =  bam::header::Header::from_template(bam_reader.header());
        let header_hashmap = bam_header.to_hashmap();
        let empty_vec = vec![];
        let read_groups = header_hashmap.get("RG").unwrap_or(&empty_vec);
        if read_groups.is_empty() {
            bail!("BAM file has no read groups (RG) tag: {}", bam_fn.to_string_lossy());
        }

        // there is only one read group
        let mut bam_sample_name: Option<String> = None;
        for read_group in read_groups.iter() {
            let rg_sample_name = match read_group.get("SM") {
                Some(s) => s,
                None => {
                    bail!("BAM file has read group with no sample name (SM) tag: {}", bam_fn.to_string_lossy());
                }
            };
            match bam_sample_name.as_ref() {
                Some(s) => {
                    if rg_sample_name != s {
                        bail!("BAM file with multiple sample reads groups detected, this is not supported: {}", bam_fn.to_string_lossy());
                    }
                },
                None => {
                    bam_sample_name = Some(rg_sample_name.clone());
                }
            };
        }

        if bam_sample_name.unwrap() == sample_name {
            sample_bams.push(bam_fn.clone());
            bam_indices.push(bam_index);
        }
    }
    Ok((sample_bams, bam_indices))
}

/// Returns true if an alignment record should be filtered out, aka ignored.
/// Main reasons for ignoring are if it is unmapped, secondary, failed QC, a duplicate, or has too low of a MAPQ.
/// # Arguments
/// * `record` - the record of interest
/// * `min_mapq` - the minimum MAPQ we allow
pub fn filter_out_alignment_record(record: &bam::Record, min_mapq: u8) -> bool {
    static FLAG_FILTER: u32 =
        htslib::BAM_FUNMAP | htslib::BAM_FSECONDARY | htslib::BAM_FQCFAIL | htslib::BAM_FDUP;

    ((record.flags() as u32) & FLAG_FILTER) != 0 || record.mapq() < min_mapq
}

/// Returns true if a VCF record can be included in phasing based on provided criteria.
/// Homozygous reference variants are always excluded unless `is_hom_allowed==true`.
/// Calls with anything missing (e.g. "./.") are also excluded.
/// Calls that match an unknown/unhandled variant type are also excluded.
/// # Arguments
/// * `record` - the variant record to check
/// * `sample_index` - the sample index, always 0 for single-sample VCFs
/// * `min_quality` - the minimum GQ for an acceptable variant
/// * `is_hom_allowed` - if true, then homozygous ALT variants are allowed, provided they meet the other criteria
/// # Errors
/// * if zygosity cannot be loaded
/// * if a call does not have a GQ tag
pub fn is_phasable_variant(record: &bcf::Record, sample_index: usize, min_quality: i32, is_hom_allowed: bool) -> Result<bool, Box<dyn std::error::Error>> {
    // check if this variant is heterozygous
    let zygosity: Zygosity = get_variant_zygosity(record, sample_index)?;
    if zygosity == Zygosity::Unknown || zygosity == Zygosity::HomozygousReference || (
        zygosity == Zygosity::HomozygousAlternate && !is_hom_allowed
    ) {
        // if unknown or homozygous reference, we definitely return false
        // if it's homozygous alternate, we also need to check if homs are allowed
        Ok(false)
    } else {
        // heterozygous, check if the variant call is of sufficient quality
        match record.format(b"GQ").integer() { // for some reason, calling .float() here will error
            Ok(all_gq) => {
                let call_quality: i32 = all_gq[sample_index][0];
                if call_quality < min_quality {
                    return Ok(false);
                }
            },
            Err(_) => {
                // usually means there is not a GQ tag, so skip this check
                // TODO: how do we long-term want to handle this for our variants?
                // trace! added mostly so clippy stops yelling at me
                trace!("Variant found without GQ tag {:?}", record);
            }
        }

        // heterozygous variant, check that the type is allowed
        let variant_type = get_variant_type(record)?;
        match variant_type {
            VariantType::Snv | 
            VariantType::Insertion |
            VariantType::Deletion |
            VariantType::Indel |
            VariantType::SvInsertion | 
            VariantType::SvDeletion |
            VariantType::TandemRepeat => { Ok(true) },

            VariantType::SvDuplication |
            VariantType::SvInversion |
            VariantType::SvBreakend |
            VariantType::Unknown=> { Ok(false) }
        }
    }
}

/// Looks at a bcf record and return the zygosity. Any "." alleles lead to Unknown zygosity results.
/// # Arguments
/// * `record` - the record to parse
/// * `sample_index` - the sample index, always 0 for single-sample VCFs
/// # Errors
/// * if rust_htslib fails to parse the genotype
/// * if the genotype field is completely empty
pub fn get_variant_zygosity(record: &bcf::Record, sample_index: usize) -> Result<Zygosity, Box<dyn std::error::Error>> {
    let all_genotypes = record.genotypes()?;
    let genotype = all_genotypes.get(sample_index);

    // if the genotype field is completely empty, something is wrong with the VCF
    if genotype.is_empty() {
        let chromosome = match record.rid() {
            Some(rid) => {
                let header = record.header();
                match header.rid2name(rid) {
                    Ok(name) => std::str::from_utf8(name).unwrap_or("FROMUTF8_ERROR"),
                    Err(_e) => "RID2NAME_ERROR"
                }
            },
            None => "NO_RID"
        };
        bail!("Encountered empty GT field for record: {}:{}", chromosome, record.pos());
    }

    let gt1 = match genotype[0] {
        GenotypeAllele::Unphased(at) => at,
        GenotypeAllele::Phased(at) => at,
        //TODO: ignore these for now, not sure how to handle it?
        GenotypeAllele::UnphasedMissing => return Ok(Zygosity::Unknown),
        GenotypeAllele::PhasedMissing => return Ok(Zygosity::Unknown)
    };

    let gt2 = if genotype.len() == 1 {
        // if the genotype has only one entry, we will just assume that gt2 is identical to gt1
        // this basically converts all single-entry genotypes into some Homozygous state
        gt1
    } else {
        match genotype[1] {
            GenotypeAllele::Unphased(at) => at,
            GenotypeAllele::Phased(at) => at,
            //TODO: ignore these for now, not sure how to handle it?
            GenotypeAllele::UnphasedMissing => return Ok(Zygosity::Unknown),
            GenotypeAllele::PhasedMissing => return Ok(Zygosity::Unknown)
        }
    };
    let zygosity = if gt1 == gt2 {
        if gt1 == 0 {
            Zygosity::HomozygousReference
        } else {
            Zygosity::HomozygousAlternate
        }
    } else {
        Zygosity::Heterozygous
    };
    Ok(zygosity)
}

/// Returns a variant type based on the alleles in the VCF.
/// # Arguments
/// * `record` - the variant record to check
pub fn get_variant_type(record: &bcf::Record) -> Result<VariantType, Box<dyn std::error::Error>> {
    // check if this has an SVTYPE field and parse into an SV type if it does
    let svtype_result = record.info("SVTYPE".as_bytes()).string();
    match svtype_result {
        Ok(svtype_option) => {
            if let Some(svtype) = svtype_option {
                // svtype is an array of strings at this point, make sure we only get one
                assert_eq!(svtype.len(), 1);
                
                // make sure these only have one ALT allele
                let num_alleles = record.alleles().len();
                assert_eq!(num_alleles, 2);
                
                let svtype_str = std::str::from_utf8(svtype[0]).unwrap();
                let sv_tag = match svtype_str {
                    "DEL" => {
                        VariantType::SvDeletion
                    },
                    "INS" => {
                        VariantType::SvInsertion
                    },
                    "DUP" => {
                        VariantType::SvDuplication
                    },
                    "INV" => {
                        VariantType::SvInversion
                    },
                    "BND" => {
                        VariantType::SvBreakend
                    },
                    _ => {
                        bail!("Unhandled SVTYPE tag: {:?}", svtype_str);
                    }
                };
                return Ok(sv_tag);
            };
        },
        Err(rust_htslib::errors::Error::BcfUndefinedTag{ tag: _ }) => {},
        Err(e) => {
            // no SVTYPE entry, so we assume it matches SNV or indel models
            bail!("Unexpected error: {:?}", e);
        }
    }

    let trid_result = record.info("TRID".as_bytes()).string();
    match trid_result {
        Ok(_trid) => {
            // we found a TRID field this is a tandem repeat
            return Ok(VariantType::TandemRepeat);
        },
        Err(rust_htslib::errors::Error::BcfUndefinedTag{ tag: _ }) => {},
        Err(e) => {
            // no SVTYPE entry, so we assume it matches SNV or indel models
            bail!("Unexpected error: {:?}", e);
        }
    }
    // TODO: we may eventually need to add a check that verifies that the only REF and ALT alleles at this point are in 
    // the normal ACGTN alphabet

    // if we have no ALT alleles and know tags to inform us, we have know idea what this is
    if record.alleles().len() <= 1 {
        return Ok(VariantType::Unknown);
    }
    
    // reference length is pulled out first, then we can look at the other alleles
    let ref_len = record.alleles()[0].len();

    // we only care about max ALT length when defining small variant type
    let max_alt_len = record.alleles().iter().skip(1)
        .map(|a| a.len())
        .max()
        .unwrap();

    Ok(if ref_len == 1 {
        if max_alt_len == 1 {
            VariantType::Snv
        } else {
            VariantType::Insertion
        }
    } else if max_alt_len == 1 {
        VariantType::Deletion
    } else {
        VariantType::Indel
    })
}

/// Defines a subset of the total reference space that is a single phasing problem or "block".
/// Each block has at least 1 read spanning from one variant to the next.
#[derive(Clone, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct PhaseBlock {
    // NOTE: order matters here because we're deriving the comparisons
    /// An index of the block, important for maintaining output order downstream.
    block_index: usize,
    /// The chromosome of the block.
    chrom: String,
    /// The chromosome index in the first VCF file
    chrom_index: u32,
    /// The coordinate of the first variant in the block, inclusive.
    start: u64,
    /// The coordinate of the last variant in the block, inclusive.
    end: u64,
    /// The total number of variants in the block so far.
    num_variants: usize,
    /// The number of variants from each VCF in this block
    vcf_index_counts: Vec<usize>,
    /// The minimum quality of variants that were included
    min_quality: i32,
    /// The sample name within the VCF that this block corresponds to
    sample_name: String,
    /// if True, then the variants in this block are meant to be skipped and left unphased
    unphased_block: bool
}

impl std::fmt::Debug for PhaseBlock {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // main purpose of custom was to munge the coordinates into a single string
        let mut result = f.debug_struct("PhaseBlock");
        result.field("block_index", &self.block_index)
            .field("coordinates", &format!("{}:{}-{}", self.chrom, self.start, self.end))
            .field("num_variants", &self.num_variants);
        if self.min_quality > 0 {
            // theres no real reason to spit this out unless it is doing something
            result.field("min_quality", &self.min_quality);
        }
        result.field("sample_name", &self.sample_name)
            .finish()
    }
}

impl PhaseBlock {
    /// Initializes a phase block with no variants
    /// # Arguments
    /// * `block_index` - the index of this block
    /// * `chrom` - the chromosome of the phase block
    /// * `chrom_index` - the chromosome index in the VCF file, for ordering
    /// * `min_quality` - the minimum quality to include a variant in this phase block
    /// * `sample_name` - the name of the sample in the VCF(s) that this block info corresponds to 
    /// * `num_vcfs` - the number of VCFs, used to generate initial empty counts Vec
    pub fn new(block_index: usize, chrom: String, chrom_index: u32, min_quality: i32, sample_name: String, num_vcfs: usize) -> PhaseBlock {
        PhaseBlock {
            block_index,
            chrom,
            chrom_index,
            start: 0,
            end: 0,
            num_variants: 0,
            vcf_index_counts: vec![0; num_vcfs],
            min_quality,
            sample_name,
            unphased_block: false
        }
    }

    pub fn get_block_index(&self) -> usize {
        self.block_index
    }

    pub fn set_block_index(&mut self, new_index: usize) {
        self.block_index = new_index;
    }

    pub fn get_chrom(&self) -> &str {
        &self.chrom
    }

    pub fn get_chrom_index(&self) -> u32 {
        self.chrom_index
    }

    pub fn get_start(&self) -> u64 {
        self.start
    }

    pub fn get_end(&self) -> u64 {
        self.end
    }

    pub fn get_num_variants(&self) -> usize {
        self.num_variants
    }

    pub fn vcf_index_counts(&self) -> &[usize] {
        &self.vcf_index_counts
    }

    pub fn get_min_quality(&self) -> i32 {
        self.min_quality
    }

    /// Returns the number of contained base pairs in the phase block.
    pub fn bp_len(&self) -> u64 {
        self.end - self.start + 1
    }

    pub fn sample_name(&self) -> &str {
        &self.sample_name
    }

    pub fn set_unphased_block(&mut self) {
        self.unphased_block = true;
    }
    
    pub fn unphased_block(&self) -> bool {
        self.unphased_block
    }

    /// Add a single-position variant to the phase block, will panic if the chromosome does not match
    /// # Arguments
    /// * `chrom` - the chromosome string
    /// * `pos` - the position of the variant
    /// * `vcf_index` - the index of the VCF this variant comes from; use 0 if only one VCF is being used
    pub fn add_locus_variant(&mut self, chrom: &str, pos: u64, vcf_index: usize) {
        assert_eq!(self.chrom, chrom, "PhaseBlock chromosomes are not equal: \"{}\" \"{}\"", self.chrom, chrom);
        //first condition unlikely to happen; second happens on an empty block
        if self.start > pos || self.num_variants == 0 {
            self.start = pos;
        }
        //most common case, extend to the right
        if self.end < pos {
            self.end = pos;
        }
        self.num_variants += 1;
        self.vcf_index_counts[vcf_index] += 1;
    }

    /// Checks if a given start/end overlaps the existing phase block
    /// # Arguments
    /// * `other_start` - the start position, inclusive
    /// * `other_end` - the end position, exclusive
    pub fn is_overlapping(&self, other_start: u64, other_end: u64) -> bool {
        let max_start = self.start.max(other_start);
        let min_end = (self.end+1).min(other_end);
        max_start < min_end 
    }
}

/// Iterator that will generate phase blocks consisting of a single "problem" to phase
pub struct PhaseBlockIterator {
    /// The index of the next block to yield
    next_block_index: usize,
    /// The primary traversal reader
    ref_vcf_readers: Vec<RefCell<bcf::IndexedReader>>,
    /// A copy of the VCF header, cached here for performance
    vcf_headers: Vec<bcf::header::HeaderView>,
    /// The name of the sample we care about in the VCF
    sample_name: String,
    /// The indices in the VCF file corresponding to `sample_name`
    sample_indices: Vec<usize>,
    /// Secondary traversals needed to figure out which variants can be phased
    bam_readers: Vec<RefCell<bam::IndexedReader>>,
    /// Index is based on the bcf::IndexedReader lookups
    chrom_index: u32,
    /// Position is as well, 0-based
    chrom_position: u64,
    /// The minimum allowed variant quality
    min_quality: i32,
    /// The minimum MAPQ to include a read
    min_mapq: u8,
    /// The minimum number of reads spanning two loci to connect them into a block
    min_spanning_reads: usize,
    /// if true, then supplemental mappings are allowed to join blocks
    allow_supplemental_joins: bool,
    /// Statistics on encountered variants while we iterate
    variant_stats: HashMap<(u32, VariantType, Zygosity), usize>
}

impl PhaseBlockIterator {
    /// Creates a new `PhaseBlockIterator` from a VCF file and collection of BAM files.
    /// # Arguments
    /// * `vcf_paths` - the VCF files to load variants from, must be zipped and indexed
    /// * `bam_paths` - the BAM files to load reads from, must be indexed
    /// * `reference_filename` - the reference genome filename
    /// * `sample_name` - the sample name in the VCF file
    /// * `min_quality` - the minimum quality to include a variant in a phase block
    /// * `min_mapq` - the minimum MAPQ to include a read
    /// * `min_spanning_reads` - the minimum number of reads that must span two adjacent variants to be joined into a phase block
    /// * `allow_supplemental_joins` - if True, supplemental mappings are used for extending blocks
    /// * `thread_pool` - a shared thread pool for BAM I/O
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        vcf_paths: &[PathBuf], bam_paths: &[PathBuf], reference_filename: &Path,
        sample_name: String,
        min_quality: i32, min_mapq: u8, min_spanning_reads: usize,
        allow_supplemental_joins: bool,
        thread_pool: &rust_htslib::tpool::ThreadPool
    ) -> Result<PhaseBlockIterator, Box<dyn std::error::Error>> {
        // needed for header() extraction
        use rust_htslib::bcf::Read;
        
        let mut ref_vcf_readers: Vec<RefCell<bcf::IndexedReader>> = vec![];
        let mut vcf_headers: Vec<bcf::header::HeaderView> = vec![];
        let mut vcf_contigs: Vec<HashSet<String>> = vec![];
        
        let mut sample_indices: Vec<usize> = vec![];

        for path in vcf_paths.iter() {
            let vcf_reader: bcf::IndexedReader = bcf::IndexedReader::from_path(path)?;
            let vcf_header: bcf::header::HeaderView = vcf_reader.header().clone();
            let ref_vcf_reader: RefCell<bcf::IndexedReader> = RefCell::new(vcf_reader);

            // first, check the sample names
            let sample_index = {
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
                        index
                    },
                    None => {
                        bail!("Sample name {:?} was not found in VCF: {:?}", sample_name, path);
                    }
                }
            };

            let contig_count: usize = vcf_header.contig_count() as usize;
            let contigs: HashSet<String> = (0..contig_count)
                .map(|i| 
                    std::str::from_utf8(
                        vcf_header.rid2name(i as u32).unwrap()
                    ).unwrap().to_string()
                )
                .collect();

            // push everything to our lists
            ref_vcf_readers.push(ref_vcf_reader);
            vcf_headers.push(vcf_header);
            vcf_contigs.push(contigs);
            sample_indices.push(sample_index);
        }

        // check that our chromosome entries are a match, if not we explode
        let first_chromosomes = &vcf_contigs[0];
        for other_chroms in vcf_contigs.iter().skip(1) {
            if first_chromosomes != other_chroms {
                bail!("Contig sets in the VCF files do not match");
            }
        }

        // open up the bam files as well
        let mut bam_readers: Vec<RefCell<bam::IndexedReader>> = vec![];
        for path in bam_paths.iter() {
            use rust_htslib::bam::Read;
            let mut bam_reader = bam::IndexedReader::from_path(path)?;
            bam_reader.set_reference(reference_filename)?;
            bam_reader.set_thread_pool(thread_pool)?;
            bam_readers.push(RefCell::new(bam_reader));
        }
        assert!(min_spanning_reads > 0);

        debug!("Sample \"{}\" VCF indices: {:?}", sample_name, sample_indices);
        
        Ok(PhaseBlockIterator {
            next_block_index: 0,
            ref_vcf_readers,
            vcf_headers,
            sample_name,
            sample_indices,
            bam_readers,
            chrom_index: 0,
            chrom_position: 0,
            min_quality,
            min_mapq,
            min_spanning_reads,
            allow_supplemental_joins,
            variant_stats: Default::default()
        })
    }

    pub fn sample_indices(&self) -> &[usize] {
        &self.sample_indices
    }
    pub fn sample_name(&self) -> &str {
        &self.sample_name
    }

    /// Retrieves variant counts from all parsed variants (whether included or not), best if used after done iterating
    pub fn variant_stats(&self) -> HashMap<(String, VariantType, Zygosity), usize> {
        let mut ret: HashMap<(String, VariantType, Zygosity), usize> = Default::default();

        for (&(chrom_index, variant_type, zygosity), &count) in self.variant_stats.iter() {
            //get the chromosome name, we iterate based on the order of the first VCF provided
            let chrom_name: String = std::str::from_utf8(
                self.vcf_headers[0].rid2name(chrom_index).unwrap()
            ).unwrap().to_string();

            ret.insert((chrom_name, variant_type, zygosity), count);
        }

        ret
    }

    /// Returns the farthest position of reads that spans chrom:pos such that _at least_ `min_read_count` reads cover that position.
    /// Currently, supplemental alignments are each handled separately.
    /// # Arguments
    /// * `chrom` - the chromosome of the locus
    /// * `pos` - the position of the locus
    fn get_longest_multispan(&self, chrom: &str, pos: u64) -> u64 {
        use bio::bio_types::genome::AbstractInterval;
        use rust_htslib::bam::Read;
        let mut span_list: Vec<u64> = vec![];
        for bam_ref in self.bam_readers.iter() {
            let mut bam = bam_ref.borrow_mut();
            bam.fetch((chrom, pos, pos+1)).unwrap();
            
            // calling .records() is what is triggering the URL warning
            for read_entry in bam.records() {
                let mut read = read_entry.unwrap();
                
                //make sure we care about the alignment
                if filter_out_alignment_record(&read, self.min_mapq) {
                    continue;
                }
                
                // see if this mapping goes farther than anything else so far
                read.cache_cigar();
                let full_range = read.range();
                // assertions always checked out, can remove
                // assert!(full_range.start == read.pos() as u64);
                // assert!(full_range.contains(&pos));
                span_list.push(full_range.end);
            }
        }

        if span_list.len() < self.min_spanning_reads {
            // this is a sentinel indicating that the range is effectively empty
            pos
        } else {
            span_list.sort();
            span_list[span_list.len() - self.min_spanning_reads]
        }
    }

    /// Returns the next position after `pos` such that _at least_ `min_read_count` reads have been found.
    /// # Arguments
    /// * `chrom` - the chromosome of the locus
    /// * `pos` - the position of the locus
    fn get_next_mapped(&self, chrom: &str, pos: u64) -> u64 {
        use rust_htslib::bam::Read;

        // TODO: This method may over-generate blocks when a user uses something other than `min_read_count = 1`.
        //       If this becomes an issue, we will likely need to rework this code a little bit to read the cigar.
        
        let mut next_positions: Vec<i64> = vec![];
        // iterate over each bam, we will cache the next `min_read_count` positions for each BAM
        for bam_ref in self.bam_readers.iter() {
            let mut bam = bam_ref.borrow_mut();
            bam.fetch((chrom, pos, i64::MAX)).unwrap();
            let mut counted: usize = 0;
            
            // calling .records() is what is triggering the URL warning
            for read_entry in bam.records() {
                let read = read_entry.unwrap();
                
                //make sure we care about the alignment
                if filter_out_alignment_record(&read, self.min_mapq) {
                    continue;
                }
                
                let start_position: i64 = read.pos();
                next_positions.push(start_position);
                counted += 1;
                if counted >= self.min_spanning_reads {
                    // we found enough reads from this BAM
                    break;
                }
            }
        }

        if next_positions.len() >= self.min_spanning_reads {
            // we found enough reads, sort them by position, then return the one at the correct location
            next_positions.sort();
            next_positions[self.min_spanning_reads-1] as u64
        } else {
            // we did not find enough reads on the rest of the chromosome, so return the max
            u64::MAX
        }

    }

    /// Returns true if there are at least `min_read_count` reads that connect from the given position back into the current phase block.
    /// # Arguments
    /// * `chrom` - the chromosome of the locus
    /// * `pos` - the position of the locus
    fn is_supplemental_overlap(&self, chrom: &str, pos: u64, phase_block: &PhaseBlock) -> bool {
        use rust_htslib::bam::Read;
        use rust_htslib::bam::record::{Aux, CigarString, Record};
        use rust_htslib::bam::record::Cigar;
        let mut overlap_count: usize = 0;
        for bam_ref in self.bam_readers.iter() {
            let mut bam = bam_ref.borrow_mut();
            bam.fetch((chrom, pos, pos+1)).unwrap();
            for read_entry in bam.records() {
                let read: Record = read_entry.unwrap();
                
                // make sure we care about the alignment
                if filter_out_alignment_record(&read, self.min_mapq) {
                    continue;
                }
                
                // check if we have any supplemental alignments
                let sa_tag: &str = match read.aux(b"SA") {
                    Ok(value) => {
                        match value {
                            Aux::String(tag) => tag,
                            _ => panic!("Unexpected tag {value:?}")
                        }
                    },
                    Err(_) => {
                        continue;
                    }
                };

                // there can be multiple, so split on the delimiter and handle each one separately
                let sa_strings: Vec<&str> = sa_tag.split_terminator(';').collect();
                for &sa_str in sa_strings.iter() {
                    // we expect exactly 6
                    let sa_frags: Vec<&str> = sa_str.split(',').collect();
                    assert_eq!(sa_frags.len(), 6);

                    let sa_chrom = sa_frags[0];
                    let sa_mapq: u8 = sa_frags[4].parse().unwrap();
                    if sa_chrom != chrom || sa_mapq < self.min_mapq {
                        // different chromosome OR mapq of the SA is too low, skip it
                        continue;
                    }
                    
                    // convert the start coordinate + CIGAR into an end coordinate
                    let sa_start: u64 = sa_frags[1].parse().unwrap();
                    let mut sa_end: u64 = sa_start;
                    let cigar: CigarString = CigarString::try_from(sa_frags[3]).unwrap();
                    for cigar_value in cigar.iter() {
                        match cigar_value {
                            Cigar::SoftClip(_) | Cigar::Ins(_) => {},
                            Cigar::Match(c_len) | 
                            Cigar::Del(c_len) | 
                            Cigar::Equal(c_len) |
                            Cigar::Diff(c_len) => {
                                sa_end += *c_len as u64;
                            },
                            _ => {
                                panic!("Unhandled cigar type: {cigar_value:?}");
                            }
                        }
                    }
                    
                    // we have the SA start and end, see if it overlaps the existing block
                    // TODO: this isn't checking for variant overlaps, can we fix that somehow? maybe instead of phase block, we pass in a variant interval tree?
                    //       this is fortunately not a major problem, any falsely connected blocks will get split on the back end
                    let overlapping: bool = phase_block.is_overlapping(sa_start, sa_end);
                    if overlapping {
                        // we found at least one overlap with the block, so increment and go to next read
                        overlap_count += 1;
                        break;
                    }
                }
            }
        }

        // return true if we got enough supplements connecting us
        overlap_count >= self.min_spanning_reads
    }
}

impl Iterator for PhaseBlockIterator {
    type Item = Result<PhaseBlock, Box<dyn std::error::Error>>;

    fn next(&mut self) -> Option<Result<PhaseBlock, Box<dyn std::error::Error>>> {
        use rust_htslib::bcf::Read;

        // make sure we still have chromosome to iterate on
        let num_contigs: u32 = self.vcf_headers[0].contig_count();
        if self.chrom_index < num_contigs {
            //get the chromosome name, we iterate based on the order of the first VCF provided
            let chrom_name: String = std::str::from_utf8(
                self.vcf_headers[0].rid2name(self.chrom_index).unwrap()
            ).unwrap().to_string();

            // initialize with an empty block containing just this chromosome
            let mut phase_block: PhaseBlock = PhaseBlock::new(
                self.next_block_index, chrom_name.clone(), self.chrom_index, self.min_quality, self.sample_name.clone(),
                self.ref_vcf_readers.len()
            );
            self.next_block_index += 1;
            
            // initalize the variant queue with one variant from each VCF, any ties in position are broken by VCF input order
            let mut variant_queue: PriorityQueue<usize, (Reverse<i64>, Reverse<usize>)> = PriorityQueue::new();
            let mut vcf_readers: Vec<_> = self.ref_vcf_readers.iter().map(|rvr| rvr.borrow_mut()).collect();
            let mut vcf_iterators: Vec<_> = vec![];
            for (vcf_index, vcf_reader) in vcf_readers.iter_mut().enumerate() {
                // fetch the corresponding chrom index for this VCF file (they are not guaranteed to match)
                let chrom_index: u32 = self.vcf_headers[vcf_index].name2rid(chrom_name.as_bytes()).unwrap();

                // fetch our position in the VCF file
                match vcf_reader.fetch(chrom_index, self.chrom_position, None) {
                    Ok(()) => {
                        // we have entries, so get the first one and queue it
                        let mut vcf_iter = vcf_reader.records().peekable();
                        let first_entry = vcf_iter.peek();
                        if let Some(record_result) = first_entry {
                            let record: &rust_htslib::bcf::Record = match record_result {
                                Ok(r) => r,
                                // we have to convert to an owned error here, and the htslib errors are not cloneable
                                Err(e) => return Some(Err(Box::new(SimpleError::from(e))))
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

            if variant_queue.is_empty() {
                // this must be an empty chromosome block because neither iterator found stuff to iterate on
                self.chrom_index += 1;
                return Some(Ok(phase_block));
            }

            let mut previous_pos: u64 = 0;
            let mut max_span: u64 = 0;
            let mut next_valid_read_pos: u64 = 0;

            while !variant_queue.is_empty() {
                // get the source of the next variant to process
                let (pop_index, pop_priority) = variant_queue.pop().unwrap();
                let sample_index = self.sample_indices[pop_index];

                // process this variant
                let record_result = vcf_iterators[pop_index].next().unwrap();
                let record = match record_result {
                    Ok(r) => r,
                    Err(e) => return Some(Err(Box::new(e)))
                };

                let variant_pos = record.pos() as u64;
                assert_eq!(variant_pos, pop_priority.0.0 as u64); // sanity check that the variant matches our position priority
                if variant_pos < self.chrom_position {
                    // this can happen when you have very very long indels that span one of our breaks
                    // we have already written though, so don't write it again
                    // skip down to variant advancement
                } else {
                    let include_variant = match is_phasable_variant(&record, sample_index, self.min_quality, false) {
                        Ok(iv) => iv,
                        Err(e) => return Some(Err(e))
                    };

                    // second condition is for variants that overlap but are before our start position
                    if include_variant {
                        //heterozygous variant found
                        if phase_block.get_num_variants() == 0 {
                            //this is a new block, first add the variant
                            phase_block.add_locus_variant(&chrom_name, variant_pos, pop_index);
                            
                            // go ahead and run the max span calculation
                            max_span = self.get_longest_multispan(&chrom_name, variant_pos);
                            if max_span == variant_pos {
                                // there are not enough reads overlapping this position, it will be unphased
                                phase_block.set_unphased_block();
                                next_valid_read_pos = self.get_next_mapped(&chrom_name, variant_pos);

                                // this just insures than any overlapping variants get in to match previous run pattern
                                max_span += 1;
                            }
                        }
                        else if max_span > variant_pos {
                            // we already found enough reads that spans _past_ this variant, so just add it
                            phase_block.add_locus_variant(&chrom_name, variant_pos, pop_index);
                        } else if phase_block.unphased_block() {
                            // if we get here, we are using logic for an unphased phase block
                            if variant_pos < next_valid_read_pos {
                                // this variant is before the next valid read, so it's automatically unphased as well
                                phase_block.add_locus_variant(&chrom_name, variant_pos, pop_index);
                            } else {
                                // this variant potentially has reads to use, so we return our unphased block
                                self.chrom_position = variant_pos;
                                return Some(Ok(phase_block));
                            }
                        } else {
                            //we check the reads from the most recent locus
                            max_span = self.get_longest_multispan(&chrom_name, previous_pos);
                            assert!(max_span != previous_pos);
                            if max_span > variant_pos {
                                //new max span connects
                                phase_block.add_locus_variant(&chrom_name, variant_pos, pop_index);
                            } else if !self.allow_supplemental_joins {
                                // no mapping spans and we are not allowing supplemental mappings to make the join
                                self.chrom_position = variant_pos;
                                return Some(Ok(phase_block));
                            } else {
                                //no *mappings* span both this new position and the most recent, check if we can find a supplemental mapping that does
                                let supplemental_overlap: bool = self.is_supplemental_overlap(&chrom_name, variant_pos, &phase_block);
                                if supplemental_overlap {
                                    // we got a supplemental mapping that works, so add this locus and go on as normal
                                    phase_block.add_locus_variant(&chrom_name, variant_pos, pop_index);
                                } else {
                                    // no overlapping mapping and no supplemental either, time to end the block
                                    self.chrom_position = variant_pos;
                                    return Some(Ok(phase_block));
                                }
                            }
                        }

                        previous_pos = variant_pos;
                    }

                    // at this point either:
                    // 1) we added the variant to the current block and are looping back around OR 
                    // 2) we did NOT add the variant to the block because it isn't phasable OR
                    // 3) we finished a block and returned out (aka, we can't get here if we just ended a block)
                    // this means that these variants are safe to add to our stats without being double counted or previously counted
                    let variant_type: VariantType = match get_variant_type(&record) {
                        Ok(vt) => vt,
                        // we have to convert to an owned error here, and the htslib errors are not cloneable
                        Err(e) => return Some(Err(e))
                    };
                    let zygosity: Zygosity = match get_variant_zygosity(&record, sample_index) {
                        Ok(z) => z,
                        // we have to convert to an owned error here, and the htslib errors are not cloneable
                        Err(e) => return Some(Err(e))
                    };

                    // update our variant stats for reporting later
                    let stats_entry = self.variant_stats.entry((self.chrom_index, variant_type, zygosity)).or_insert(0);
                    *stats_entry += 1;
                }

                // requeue from the one we popped from
                let next_entry = vcf_iterators[pop_index].peek();
                if let Some(record_result) = next_entry {
                    let record: &rust_htslib::bcf::Record = match record_result {
                        Ok(r) => r,
                        // we have to convert to an owned error here, and the htslib errors are not cloneable
                        Err(e) => return Some(Err(Box::new(SimpleError::from(e))))
                    };
                    let position: i64 = record.pos();
                    variant_queue.push(pop_index, (Reverse(position), Reverse(pop_index)));
                };
            }

            //we have reached the end of the current chromosome, reset to next chromosome and return what we have
            self.chrom_index += 1;
            self.chrom_position = 0;
            Some(Ok(phase_block))
        } else {
            // no chromosomes left to iterate on
            None
        }
    }
}

/// Iterator over multiple phase blocks iterators.
/// Output blocks are ordered by (chromosome, start_position, end_position) and re-numbered to reflect traversal order.
pub struct MultiPhaseBlockIterator {
    /// The internal iterators we use
    sub_iterators: Vec<PhaseBlockIterator>,
    /// The priority queue for the phase blocks at the front of each iterator
    phase_block_queue: PriorityQueue<(usize, PhaseBlock), PhaseBlockPriority>,
    /// The combined block index
    joint_block_index: usize,
}

impl MultiPhaseBlockIterator {
    /// Creates a new iterator from a vector of sub-iterators, each one tied to a specific sample.
    /// # Arguments
    /// * `sub_iterators` - the original PhaseBlockIterators that this will wrap
    /// # Errors
    /// * if any sub-iterators generate errors while iterating
    pub fn new(mut sub_iterators: Vec<PhaseBlockIterator>) -> Result<MultiPhaseBlockIterator, Box<dyn std::error::Error>> {
        let mut phase_block_queue: PriorityQueue<(usize, PhaseBlock), PhaseBlockPriority> = PriorityQueue::new();

        for (index, iterator) in sub_iterators.iter_mut().enumerate() {
            let next_value = iterator.next();
            match next_value {
                Some(result) => {
                    let first_block: PhaseBlock = result?;
                    let block_priority = Self::get_block_priority(&first_block);
                    phase_block_queue.push((index, first_block), block_priority);
                },
                None => {
                    // first block is empty, which is weird but technically allowed
                    warn!("First block in iterator {} was empty.", index);
                }
            };
        }

        Ok(MultiPhaseBlockIterator { 
            sub_iterators, 
            phase_block_queue,
            joint_block_index: 0
        })
    }

    /// Retrieves variant counts from all parsed variants (whether included or not) across all samples.
    /// This is best if done when iteration is finished.
    /// Returns a hashmap where key is (sample_name, chromosome, variant_type, zygosity) and value is a count.
    pub fn variant_stats(&self) -> HashMap<(String, String, VariantType, Zygosity), usize> {
        // key = (sample_name, chromosome, variant_type, zygosity); value = count
        let mut ret: HashMap<(String, String, VariantType, Zygosity), usize> = Default::default();
        for pbi in self.sub_iterators.iter() {
            let sample_name = pbi.sample_name().to_string();
            let pbi_stats = pbi.variant_stats();
            for ((chrom, vt, zyg), count) in pbi_stats.into_iter() {
                ret.insert((sample_name.clone(), chrom, vt, zyg), count);
            }
        }
        ret
    }

    /// Returns the block priority for a phase block
    /// # Arguments
    /// * `phase_block` - the block to calculate priority for
    fn get_block_priority(phase_block: &PhaseBlock) -> PhaseBlockPriority {
        Reverse(
            (
                phase_block.get_chrom_index(),
                phase_block.get_start(),
                phase_block.get_end()
            )
        )
    }
}

impl Iterator for MultiPhaseBlockIterator {
    type Item = Result<PhaseBlock, Box<dyn std::error::Error>>;

    fn next(&mut self) -> Option<Result<PhaseBlock, Box<dyn std::error::Error>>> {
        let pq_next = self.phase_block_queue.pop();
        match pq_next {
            Some(((source_index, mut phase_block), _priority)) => {
                // get the next block and put it on the queue
                let next_item = self.sub_iterators[source_index].next();
                if let Some(next_result) = next_item {
                    // we have more in this queue, add it to the priority queue
                    let next_block = match next_result {
                        Ok(b) => b,
                        Err(e) => {
                            // sub-queue error, propagate it up the chain
                            return Some(Err(e));
                        }
                    };
                    let next_priority = Self::get_block_priority(&next_block);
                    self.phase_block_queue.push((source_index, next_block), next_priority);
                };

                // we need to update the block index based on the joint values
                phase_block.set_block_index(self.joint_block_index);
                self.joint_block_index += 1;

                // finally send back the block
                Some(Ok(phase_block))
            },
            None => {
                None
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // fortunately, none of the BAM checks actually look at the provided reference file contents 
    const LOCAL_REFERENCE: &str = "./test_data/test_reference.fa";
    
    #[test]
    fn test_get_vcf_samples() {
        let header_only_vcf: PathBuf = "./test_data/header_only.vcf.gz".into();
        let expected_samples: Vec<String> = vec![
            "HG001".to_string(), "HG002_30x".to_string(), "HG005_30x".to_string()
        ];
        let samples = get_vcf_samples(&header_only_vcf).unwrap();
        assert_eq!(expected_samples, samples);
    }

    #[test]
    fn test_get_sample_bams() {
        let sample_name: String = "HG002-rep1".to_string();
        let all_bams: Vec<PathBuf> = vec![
            "./test_data/header_only.bam".into(),
            "./test_data/multi_smrtcell.bam".into()
        ];
        
        let (bams_found, bam_indices) = get_sample_bams(
            &all_bams,
            &sample_name,
            &PathBuf::from(LOCAL_REFERENCE)
        ).unwrap();
        assert_eq!(all_bams, bams_found);
        assert_eq!(vec![0, 1], bam_indices);
    }

    #[test]
    fn test_multisample_bam() {
        let sample_name: String = "HG002-rep1".to_string();
        let all_bams: Vec<PathBuf> = vec![
            "./test_data/multisample.bam".into()
        ];
        let result = get_sample_bams(
            &all_bams,
            &sample_name,
            &PathBuf::from(LOCAL_REFERENCE)
        );

        // should have an error with the following message
        assert!(result.is_err());
        let expected_error_string = "BAM file with multiple sample reads groups detected, this is not supported: ./test_data/multisample.bam".to_string();
        assert_eq!(expected_error_string, result.err().unwrap().to_string());
    }
}