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
