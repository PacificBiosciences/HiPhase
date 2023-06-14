# v0.10.1
## Fixed
- Corrected a panic caused by reference mismatches to produce a human-readable error message (Resolves #12)
- Fixed an issue parsing unexpected CIGAR strings for local re-alignment and SV deletions
- Added a check for gzip encoding and tabix index for input VCF file to prevent uninterpretable segfault in htslib (Resolves #13)

# v0.10.0
## Changes
- Adds support for multi-threaded BAM I/O via shared thread pools. Internal testing with default parameters showed ~40% reduction in wall-clock time when haplotagging is enabled with 16 threads. User experience will vary based on file type, disk type, and I/O contention. Resolves #9.
- Adds `--io-threads` option for greater control over the number of threads used in the thread pool. Defaults to the same number provided to `--threads`.

# v0.9.0
## Changes
- Adds support for CRAM input and output files
- The reference file (`--reference`) is now a required parameter for HiPhase. This prevents users from accidentally forgetting the reference file which leads to drastically different performance.

## Fixed
- Corrected an error where bgzipped reference files were not correctly loading

# v0.8.1
## Changes
- Adds `--ignore-read-groups` flag that will skip the read group checks for provided BAM files. This option can only be used to phase one dataset at a time. Resolves #3.
- Adds `--haplotag-file <FILE>` option that will create a TSV/CSV file containing read haplotag information. Resolves #4.

## Fixed
- Corrected some typos in CLI help menu
- Changes HP tag in output BAM files from an i32 to a u8

# v0.8.0
## Changes
Multi-sample VCF inputs are now supported: 
- Users can provide VCFs with multiple samples; all samples being phased must appear in all provided VCF files
- Each sample will be phased independently of the other samples, and results placed in a phased, multi-sample VCF
- Only single sample BAM files are supported, and the sample name *must* be present in the read group (RG) sample (SM) tag and match the VCF sample identifier
- Haplotagging with multi-sample VCFs is also supported

# v0.7.3
## Fixed
- Adjust linux release build process to improve portability.  This release should be functionally identical to v0.7.2 so no update is needed for users without binary issues.
