
/// Contains writer for phase block stats, both the phase blocks themselves and the summary
pub mod block_stats;
/// Contains the writer for haplotag results
pub mod haplotag_writer;
/// Contains writer for BAM files
pub mod ordered_bam_writer;
/// Contains writer for VCF files
pub mod ordered_vcf_writer;
/// Contains writer for phasing statistics for underlying algorithms
pub mod phase_stats;
/// Contains additional VCF utilities
pub mod vcf_util;