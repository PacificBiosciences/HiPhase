
use clap::Parser;
use chrono::Datelike;
use flate2::bufread::MultiGzDecoder;
use lazy_static::lazy_static;
use log::{error, info, warn};
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

lazy_static! {
    /// Stores the full version string we plan to use.
    /// # Examples
    /// * `0.11.0-6bb9635-dirty` - while on a dirty branch
    /// * `0.11.0-6bb9635` - with a fresh commit
    pub static ref FULL_VERSION: String = format!("{}-{}", env!("CARGO_PKG_VERSION"), env!("VERGEN_GIT_DESCRIBE"));
}

#[derive(Clone, Parser)]
#[clap(author, 
    version = &**FULL_VERSION,
    about, 
    after_help = format!("Copyright (C) 2004-{}     Pacific Biosciences of California, Inc.
This program comes with ABSOLUTELY NO WARRANTY; it is intended for
Research Use Only and not for use in diagnostic procedures.", chrono::Utc::now().year()))]
pub struct Settings {
    /// Input alignment file in BAM format.
    #[clap(required = true)]
    #[clap(short = 'b')]
    #[clap(long = "bam")]
    #[clap(value_name = "BAM")]
    #[clap(help_heading = Some("Input/Output"))]
    pub bam_filenames: Vec<PathBuf>,

    /// Output haplotagged alignment file in BAM format.
    #[clap(short = 'p')]
    #[clap(long = "output-bam")]
    #[clap(value_name = "BAM")]
    #[clap(help_heading = Some("Input/Output"))]
    pub output_bam_filenames: Vec<PathBuf>,

    /// Input variant file in VCF format.
    #[clap(required = true)]
    #[clap(short = 'c')]
    #[clap(long = "vcf")]
    #[clap(value_name = "VCF")]
    #[clap(help_heading = Some("Input/Output"))]
    pub vcf_filenames: Vec<PathBuf>,

    /// Output phased variant file in VCF format.
    #[clap(required = true)]
    #[clap(short = 'o')]
    #[clap(long = "output-vcf")]
    #[clap(value_name = "VCF")]
    #[clap(help_heading = Some("Input/Output"))]
    pub output_vcf_filenames: Vec<PathBuf>,

    /// Reference FASTA file
    #[clap(required = true)]
    #[clap(short = 'r')]
    #[clap(long = "reference")]
    #[clap(value_name = "FASTA")]
    #[clap(help_heading = Some("Input/Output"))]
    pub reference_filename: PathBuf,

    /// Sample name to phase within the VCF (default: first sample)
    #[clap(short = 's')]
    #[clap(long = "sample-name")]
    #[clap(value_name = "SAMPLE")]
    #[clap(help_heading = Some("Input/Output"))]
    pub sample_names: Vec<String>,

    /// Ignore BAM file read group IDs
    #[clap(long = "ignore-read-groups")]
    #[clap(help_heading = Some("Input/Output"))]
    pub ignore_read_groups: bool,

    /// Output summary phasing statistics file (optional, csv/tsv)
    #[clap(long = "summary-file")]
    #[clap(value_name = "FILE")]
    #[clap(help_heading = Some("Input/Output"))]
    pub summary_filename: Option<PathBuf>,

    /// Output algorithmic statistics file (optional, csv/tsv)
    #[clap(long = "stats-file")]
    #[clap(value_name = "FILE")]
    #[clap(help_heading = Some("Input/Output"))]
    pub stats_filename: Option<PathBuf>,

    /// Output blocks file (optional, csv/tsv)
    #[clap(long = "blocks-file")]
    #[clap(value_name = "FILE")]
    #[clap(help_heading = Some("Input/Output"))]
    pub blocks_filename: Option<PathBuf>,

    /// Output haplotag file (optional, csv/tsv)
    #[clap(long = "haplotag-file")]
    #[clap(value_name = "FILE")]
    #[clap(help_heading = Some("Input/Output"))]
    pub haplotag_filename: Option<PathBuf>,

    /// Number of threads for BAM I/O (default: minimum of `--threads` or `4`)
    #[clap(long = "io-threads")]
    #[clap(value_name = "THREADS")]
    #[clap(help_heading = Some("Input/Output"))]
    pub io_threads: Option<usize>,

    /// Output .csi indices instead of .tbi/.bai
    #[clap(long = "csi-index")]
    #[clap(help_heading = Some("Input/Output"))]
    pub csi_index: bool,

    /// Number of threads to use for phasing.
    #[clap(short = 't')]
    #[clap(long = "threads")]
    #[clap(value_name = "THREADS")]
    #[clap(default_value = "1")]
    pub threads: usize,

    /// Enable verbose output.
    #[clap(short = 'v')]
    #[clap(long = "verbose")]
    #[clap(action = clap::ArgAction::Count)]
    pub verbosity: u8,

    /// Sets a minimum genotype quality (GQ) value to include a variant in the phasing
    #[clap(long = "min-vcf-qual")]
    #[clap(value_name = "GQ")]
    #[clap(default_value = "0")]
    #[clap(help_heading = Some("Variant Filtering"))]
    pub min_variant_quality: i32,

    /// Sets a minimum MAPQ to include a read in the phasing
    #[clap(long = "min-mapq")]
    #[clap(value_name = "MAPQ")]
    #[clap(default_value = "5")]
    #[clap(help_heading = Some("Mapping Filtering"))]
    pub min_mapping_quality: u8,

    /// Sets a minimum number of matched variants required for a read to get included in the scoring
    #[clap(long = "min-matched-alleles")]
    #[clap(value_name = "COUNT")]
    #[clap(default_value = "2")]
    #[clap(help_heading = Some("Mapping Filtering"))]
    pub min_matched_alleles: usize,

    /// Sets a minimum number of reads to span two adjacent variants to join a phase block
    #[clap(long = "min-spanning-reads")]
    #[clap(value_name = "READS")]
    #[clap(default_value = "1")]
    #[clap(help_heading = Some("Phase Block Generation"))]
    pub min_spanning_reads: usize,

    /// Disables the use of supplemental mappings to join phase blocks
    #[clap(long = "no-supplemental-joins")]
    #[clap(help_heading = Some("Phase Block Generation"))]
    pub disable_supplemental_joins: bool,

    /// Enables the phasing and haplotagging of singleton phase blocks
    #[clap(long = "phase-singletons")]
    #[clap(help_heading = Some("Phase Block Generation"))]
    pub phase_singletons: bool,

    /// Sets a maximum reference buffer for local realignment
    #[clap(long = "max-reference-buffer")]
    #[clap(value_name = "LENGTH")]
    #[clap(default_value = "15")]
    #[clap(help_heading = Some("Allele Assignment"))]
    pub reference_buffer: usize,

    /// Enables global realignment with a maximum allowed CPU time before fallback to local realignment
    #[clap(long = "global-realignment-cputime")]
    #[clap(value_name = "SECONDS")]
    #[clap(default_value = "0.0")]
    #[clap(help_heading = Some("Allele Assignment"))]
    pub global_realign_cputime: f32,

    /// Sets a pruning threshold on global realignment, set to 0 to disable pruning
    #[clap(long = "global-pruning-distance")]
    #[clap(value_name = "LENGTH")]
    #[clap(default_value = "500")]
    #[clap(help_heading = Some("Allele Assignment"))]
    pub wfa_prune_distance: usize,

    /// Sets the minimum queue size for the phasing algorithm
    #[clap(long = "phase-min-queue-size")]
    #[clap(value_name = "SIZE")]
    #[clap(default_value = "1000")]
    #[clap(help_heading = Some("Phasing"))]
    pub phase_min_queue_size: usize,
    
    /// Sets the queue size increment per variant in a phase block
    #[clap(long = "phase-queue-increment")]
    #[clap(value_name = "SIZE")]
    #[clap(default_value = "3")]
    #[clap(help_heading = Some("Phasing"))]
    pub phase_queue_increment: usize,

    /// Skips a number of blocks (debug only); non-0 values will cause an error on VCF output
    #[clap(long = "skip")]
    #[clap(hide = true)]
    #[clap(default_value = "0")]
    pub skip_blocks: usize,

    /// Take a number of blocks (debug only); non-0 values will cause an error on VCF output
    #[clap(long = "take")]
    #[clap(hide = true)]
    #[clap(default_value = "0")]
    pub take_blocks: usize,
}

/// Checks if a file exists and will otherwise exit
/// # Arguments
/// * `filename` - the file path to check for
/// * `label` - the label to use for error messages
fn check_required_filename(filename: &Path, label: &str) {
    if !filename.exists() {
        error!("{} does not exist: \"{}\"", label, filename.display());
        std::process::exit(exitcode::NOINPUT);
    } else {
        info!("{}: \"{}\"", label, filename.display());
    }
}

/// Checks if the VCF file exists, is bgzipped, and has an index. If it fails any of those, this will exit.
/// # Argument
/// * `filename` - the VCF file path to check
/// * `label` - the label to use for error messages
fn check_required_vcf(filename: &Path, label: &str) {
    // first check the filename normally
    check_required_filename(filename, label);

    // now we need to check that this is a bgzipped file by just trying to read a little bit of it
    // NOTE: if the user generates a gzip file (as opposed to bgzip), this will still pass :(
    //       in theory, indexing checks should fail
    let vcf_file: File = File::open(filename).unwrap();
    let file_reader = BufReader::new(vcf_file);
    let mut gz_decoder = MultiGzDecoder::new(file_reader);
    let mut small_buffer: [u8; 10] = [0; 10];
    match gz_decoder.read(&mut small_buffer) {
        Ok(_) => {},
        Err(e) => {
            if e.to_string() == "invalid gzip header" {
                error!("Error while checking {filename:?}: {e}; is the VCF bgzipped?");
            } else {
                error!("Error while checking {filename:?}: {e}");
            }
            std::process::exit(exitcode::IOERR);
        }  
    };

    // finally, verify that an index file exists, should just be tbi and csi
    let known_indices = ["tbi", "csi"];
    let mut index_found: bool = false;
    for &ki in known_indices.iter() {
        let mut extension_path = filename.to_owned()
            .into_os_string();
        extension_path.push(format!(".{ki}"));
        let extension_path: PathBuf = PathBuf::from(extension_path);
        index_found |= extension_path.exists();
    }
    if !index_found {
        error!("Error while checking {filename:?}: no tabix index found (.tbi or .csi)");
        std::process::exit(exitcode::NOINPUT);
    }

}

pub fn get_raw_settings() -> Settings {
    Settings::parse()
}

/// Do some additional checks here, we may increase these as we go.
/// Also can modify settings if needed since we're passing it around.
/// # Arguments
/// * `settings` - the raw settings, nothing has been checked other than what clap does for us.
pub fn check_settings(mut settings: Settings) -> Settings {
    //check for any of our required files
    for filename in settings.bam_filenames.iter() {
        check_required_filename(filename, "Alignment file");
    }
    for filename in settings.vcf_filenames.iter() {
        check_required_vcf(filename, "Variant file");
    }

    // make sure the number of inputs and outputs are identical
    if settings.vcf_filenames.len() != settings.output_vcf_filenames.len() {
        error!("Detected {} input VCFs and {} output VCFs, these must be equal", settings.vcf_filenames.len(), settings.output_vcf_filenames.len());
        std::process::exit(exitcode::USAGE);
    }

    // if we have any phased BAM outputs, make sure we have one for each file
    if !settings.output_bam_filenames.is_empty() && settings.bam_filenames.len() != settings.output_bam_filenames.len() {
        error!("Detected {} input BAMs and {} output BAMs, these must be equal", settings.bam_filenames.len(), settings.output_bam_filenames.len());
        std::process::exit(exitcode::USAGE);
    }

    // check optional files
    check_required_filename(&settings.reference_filename, "Reference file");

    // 0 is just a sentinel for everything
    if settings.take_blocks == 0 {
        settings.take_blocks = usize::MAX;
    }
    if settings.wfa_prune_distance == 0 {
        settings.wfa_prune_distance = usize::MAX;
    }

    // 0 doesn't make sense, so lets just error proof it up to 1
    if settings.min_spanning_reads == 0 {
        settings.min_spanning_reads = 1;
    }
    if settings.min_matched_alleles == 0 {
        settings.min_matched_alleles = 1;
    }

    // if this is not specified, then set it to the same as processing
    if settings.io_threads.is_none() {
        // setting to the same as threads generates some issues with no real benefit
        // 4 is a happy default, users can override if needed
        settings.io_threads = Some(settings.threads.min(4));
    }

    // dump stuff to the logger
    info!("Minimum call quality: {}", settings.min_variant_quality);
    info!("Minimum mapping quality: {}", settings.min_mapping_quality);
    info!("Minimum matched alleles: {}", settings.min_matched_alleles);
    if settings.min_matched_alleles > 2 {
        warn!("Setting the minimum matched alleles > 2 has not been tested.")
    }
    info!("Minimum spanning reads: {}", settings.min_spanning_reads);
    info!("Supplemental mapping block joins: {}", if settings.disable_supplemental_joins { "DISABLED" } else { "ENABLED" });
    info!("Phase singleton blocks: {}", if settings.phase_singletons { "ENABLED" } else { "DISABLED" });
    info!("Local re-alignment maximum reference buffer: +-{} bp", settings.reference_buffer);
    if settings.global_realign_cputime == 0.0 {
        info!("Global re-alignment: DISABLED");
    } else {
        info!("Global re-alignment CPU time: {} seconds", settings.global_realign_cputime);
        if settings.wfa_prune_distance == usize::MAX {
            info!("Global prune distance: DISABLED");
        } else {
            info!("Global prune distance: {}", settings.wfa_prune_distance);
        }
    }
    info!("Processing threads: {}", settings.threads);
    info!("I/O threads: {}", settings.io_threads.unwrap());
    if settings.csi_index {
        info!("CSI indexing: enabled");
    }

    //send the settings back
    settings
}
