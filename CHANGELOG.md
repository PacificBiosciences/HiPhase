# v1.4.5
## Fixed
* Fixed an error where BAM phase tags were not always properly removed prior to re-tagging, leading to a run-time error and exit

# v1.4.4
## Fixed
* Fixed an error where phasing information that was present in input files would be copied through to output files if it was not overwritten by HiPhase phasing results. HiPhase will now automatically remove this phasing information to prevent accidental mixing of phase results. 
  * For VCF files, any unphased genotypes will be switched to unphased and sorted by allele index (e.g. 1|0 -> 0/1). The "FORMAT:PS" and "FORMAT:PF" tags will either be removed entirely if the whole record is unphased or set to "." for partially phased records.
  * For BAM files, the "HP" and "PS" tags will be removed for any unphased records.

# v1.4.3
## Fixed
* Replaced a panic caused by a chromosome appearing in a VCF but not in the BAM file with a more descriptive error message
* Fixed an error caused by a multi-sample VCF with a mixture of haploid and diploid genotypes

# v1.4.2
## Changes
* Removes a 1 basepair shift from tandem repeat region calculation to support anchor base changes in TRGT v1.0.0; internal results are nearly identical before and after this change

# v1.4.1
## Changes
* Reclassifies warnings during VCF writing to debug
* Adds a section to the quickstart guide on resource requirements

# v1.4.0
## Changes
* **Major changes to dual-mode allele assignment:** Prior to this version, global realignment would revert to local realignment if the CPU cost (in seconds) exceeded a user provided threshold. While this was useful for fast-tracking noisy phase blocks, it could lead to non-deterministic output as CPU costs can vary. The thresholding has been reworked such that global realignment will revert to local realignment *for an individual mapping* if the edit distance exceeds a user provided threshold (default: 500). Additionally, global realignment will revert to local realignment *for the remainder of a putative phase block* if too many reads have reverted to local realignment (default: 50%, minimum number of failures: 50 mappings). This has the following downstream impact on results:
  * All results from HiPhase are **fully deterministic** from run to run.
  * Baseline quality scores for local realignment have been adjusted to scale at the same relative ratios as those from global realignment. 
    * When running HiPhase on _only_ small variants (e.g., local realignment mode only), this tended to slightly increase the number of switch flip errors relative to v1.3.0.
    * When running HiPhase on small, structural, and tandem repeat variants (recommended), we observed a small decrease in switch flip errors relative to v1.3.0.
  * Relative to v1.3.0, we observed reduced run-time costs for all tests (~25% reduction in both CPU time and wall-clock time, on average).
  * The number of mappings processed through global/local realignment are now tracked in the `--stats-file`.
* **Global realignment is now on by default**, reflecting our overall recommended usage of HiPhase. This can be disabled with the `--disable-global-realignment` option.
* **CLI changes:** The CLI has been updated to reflect the above algorithmic changes. These new CLI options have been added to reflect the changes:
  * `--disable-global-realignment` - This option will disable all global realignments; it is recommended if only small variant files are available for phasing
  * `--global-realignment-max-ed <DISTANCE>` - Controls the maximum allowed edit distance before reverting an individual mapping to local realignment (default: 500)
  * `--max-global-failure-ratio <FRAC>` - Controls the maximum allowed failure rates for global realignment before reverting the rest of the phase block to local realignment (default: 50%)
  * `--global-failure-count <COUNT>` - Controls the minimum number of failures required before the failure rate check is enabled (default: 50)
  * `--global-realignment-cputime <SECONDS>` - **Deprecated**, this option is now hidden on the CLI. It will produce a warning if used but has no impact on the downstream results.

# v1.3.0
## Changes
- Relaxes the requirements for SV deletion and insertion events such that they no longer require an alternate or reference allele, respectively, to have length 1

## Internal changes
- The interface for variant creation was modified to reduce panics from invalid variant construction. This modification changes all the return types for the various `Variant::new*(...)` functions from `Variant` to `Result<Variant, VariantError>`.

## Fixed
- SV events with a placeholder ALT sequence (e.g., \<DEL\>, \<INS\>) are now properly ignored by HiPhase instead of creating an error.

# v1.2.1
## Fixed
- Fixed [a rare issue](https://github.com/PacificBiosciences/pbbioconda/issues/640) where reference alleles with stripped IUPAC codes were throwing errors due to reference mismatch
- Fixed an issue where variants preceding a GraphWFA region were not ignored, potentially leading to aberrant graph structure

# v1.2.0
## Changes
- Added an option (`--csi-index`) to output .csi index files instead of .tbi/.bai
- HiPhase will now generate an error if it detects a large reference chromosome and .csi indexing is not enabled

# v1.1.0
## Changes
- Updated the way unphaseable variants are placed into artificial phase blocks. Instead of unphased singletons, these are now grouped into larger unphased blocks. Phased variants are unaffected by this change.
  - This has a major impact on I/O churning, especially in centromeric and telomeric regions. Internal tests show up to 60% reduction in wall clock time with 16 threads accompanied by up to 30% reduction in CPU time depending on the phasing mode.
  - There are now fewer total phase blocks reported due to the singleton collapse. This primarily impacts the statistic tracking optional outputs of HiPhase.

## Fixed
- Fixed a related issue where unphased singleton blocks were falsely tagged with `TR_OVERLAP`
- Fixed an issue where the default `--io-threads` frequently caused issues on long-running jobs. This now defaults to `min(--threads, 4)`, but can still be specified with `--io-threads`.

# v1.0.0
## Changes
- Added support for tandem repeat calls from [TRGT](https://github.com/PacificBiosciences/trgt); minimum supported version - v0.5.0
- VCF index files (.tbi) are now automatically generated by HiPhase
- The reference FASTA file is now a **required** parameter for HiPhase; this change prevents sub-optimal performance caused by forgetting to specify the reference genome
- HiPhase source code is now released under new LICENSE file

## Fixed
- Fixed an infinite looping error caused by variants with the same start POS but no overlapping mappings

# v0.10.2
## Fixed
- Replaced a panic caused by large WFA edit distance with a human-readable error message (Resolves #15)

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
