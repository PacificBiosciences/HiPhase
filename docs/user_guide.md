# User Guide
Table of contents:

* [Quickstart](#quickstart)
* [Common uses cases](#common-use-cases)
* [Supported upstream processes](#supported-upstream-processes)
* [Output files](#output-files)

# Quickstart
```bash
hiphase \
    --bam {IN_BAM} \
    --vcf {IN_VCF} \
    --output-vcf {OUT_VCF} \
    --reference {REFERENCE} \
    --threads {THREADS}
```

Parameters:
* `--bam {IN_BAM}` - path to a BAM file containing reads only from the sample that is being phased, this option can be specified multiple times
* `--vcf {IN_VCF}` - path to a VCF file containing the variants to phase, this option can be specified multiple times
* `--output-vcf {OUT_VCF}` - path to the output VCF that will contain the phased variants, this option must be specified the same number of times as `--vcf` 
* `--reference {REFERENCE}` - a FASTA file containing the reference genome, gzip allowed; the reference genome file is optional, but _very strongly_ recommended for optimal performance
* `--threads {THREADS}` - number of threads to use for phasing (default: 1)

## Quickstart Example
This example also includes some auxiliary outputs that may be useful for tracking phase result statistics.
For details on each, refer to the [Output files](#output-files) section.

```
hiphase \
    --threads 16 \
    --reference human_GRCh38_no_alt_analysis_set.fasta \
    --bam m64109_200805_204709.GRCh38.bam \
    --bam m64109_200807_075817.GRCh38.bam \
    --bam m64109_200808_191025.GRCh38.bam \
    --bam m64109_200810_062248.GRCh38.bam \
    --bam m64109_200813_162416.GRCh38.bam \
    --bam m64109_200815_033514.GRCh38.bam \
    --vcf HG001.GRCh38.deepvariant.vcf.gz \
    --output-vcf HG001.GRCh38.deepvariant.phased.vcf.gz \
    --stats-file HG001.stats.csv \
    --blocks-file HG001.blocks.tsv \
    --summary-file HG001.summary.tsv
        
[2023-02-02T20:13:06.152Z INFO  hiphase::cli] Alignment file: "m64109_200805_204709.GRCh38.bam"
[2023-02-02T20:13:06.152Z INFO  hiphase::cli] Alignment file: "m64109_200807_075817.GRCh38.bam"
[2023-02-02T20:13:06.152Z INFO  hiphase::cli] Alignment file: "m64109_200808_191025.GRCh38.bam"
[2023-02-02T20:13:06.152Z INFO  hiphase::cli] Alignment file: "m64109_200810_062248.GRCh38.bam"
[2023-02-02T20:13:06.152Z INFO  hiphase::cli] Alignment file: "m64109_200813_162416.GRCh38.bam"
[2023-02-02T20:13:06.152Z INFO  hiphase::cli] Alignment file: "m64109_200815_033514.GRCh38.bam"
[2023-02-02T20:13:06.153Z INFO  hiphase::cli] Variant file: "HG001.GRCh38.deepvariant.vcf.gz"
[2023-02-02T20:13:06.153Z INFO  hiphase::cli] Reference file: "human_GRCh38_no_alt_analysis_set.fasta"
[2023-02-02T20:13:06.153Z INFO  hiphase::cli] Minimum call quality: 0
[2023-02-02T20:13:06.153Z INFO  hiphase::cli] Minimum mapping quality: 0
[2023-02-02T20:13:06.153Z INFO  hiphase::cli] Minimum matched alleles: 2
[2023-02-02T20:13:06.153Z INFO  hiphase::cli] Minimum spanning reads: 1
[2023-02-02T20:13:06.153Z INFO  hiphase::cli] Supplemental mapping block joins: ENABLED
[2023-02-02T20:13:06.153Z INFO  hiphase::cli] Local re-alignment maximum reference buffer: +-15 bp
[2023-02-02T20:13:06.153Z INFO  hiphase::cli] Global re-alignment: DISABLED
[2023-02-02T20:13:06.153Z INFO  hiphase::data_types::reference_genome] Loading "human_GRCh38_no_alt_analysis_set.fasta"...
[2023-02-02T20:13:26.049Z INFO  hiphase::data_types::reference_genome] Finished loading 195 contigs.
[2023-02-02T20:13:26.338Z INFO  hiphase] Starting job pool with 16 threads...
[2023-02-02T20:13:38.683Z INFO  hiphase] Generated 100 phase blocks, latest block: PhaseBlock { block_index: 99, coordinates: "chr1:21871257-21871257", num_variants: 1, sample_name: "HG001" }
[2023-02-02T20:13:50.615Z INFO  hiphase] Generated 200 phase blocks, latest block: PhaseBlock { block_index: 199, coordinates: "chr1:45726723-45769054", num_variants: 8, sample_name: "HG001" }
[2023-02-02T20:13:58.233Z INFO  hiphase] Generated 300 phase blocks, latest block: PhaseBlock { block_index: 299, coordinates: "chr1:62628086-62628086", num_variants: 1, sample_name: "HG001" }
[2023-02-02T20:14:09.075Z INFO  hiphase] Generated 400 phase blocks, latest block: PhaseBlock { block_index: 399, coordinates: "chr1:87116613-87231282", num_variants: 135, sample_name: "HG001" }
[2023-02-02T20:14:15.183Z INFO  hiphase] Generated 500 phase blocks, latest block: PhaseBlock { block_index: 499, coordinates: "chr1:100648818-100648818", num_variants: 1, sample_name: "HG001" }
[2023-02-02T20:14:44.050Z INFO  hiphase] Generated 600 phase blocks, latest block: PhaseBlock { block_index: 599, coordinates: "chr1:143321568-143350559", num_variants: 22, sample_name: "HG001" }
[2023-02-02T20:14:54.843Z INFO  hiphase] Generated 700 phase blocks, latest block: PhaseBlock { block_index: 699, coordinates: "chr1:164168479-164453849", num_variants: 295, sample_name: "HG001" }
[2023-02-02T20:15:01.180Z INFO  hiphase] Received results for 100 phase blocks: 1.0544 blocks/sec, 228.1697 hets/sec, writer waiting on block 78
[2023-02-02T20:15:06.073Z INFO  hiphase] Generated 800 phase blocks, latest block: PhaseBlock { block_index: 799, coordinates: "chr1:188271650-188271650", num_variants: 1, sample_name: "HG001" }
[2023-02-02T20:15:13.030Z INFO  hiphase] Received results for 200 phase blocks: 1.8745 blocks/sec, 419.1398 hets/sec, writer waiting on block 178
...
[2023-02-02T20:38:21.104Z INFO  hiphase] Received results for 12673 / 12674 phase blocks: 8.4782 blocks/sec, 2067.6570 hets/sec, writer waiting on block 9876
...
[2023-02-02T20:40:11.266Z INFO  hiphase] Saving all blocks to "HG001.blocks.tsv"...
[2023-02-02T20:40:11.340Z INFO  hiphase] Saving summary block statistics to "HG001.summary.tsv"...
[2023-02-02T20:40:11.347Z INFO  hiphase] All phase blocks finished successfully after 1605.008968506 seconds.
```

# Common use cases
## Joint phasing small variants and SVs
To *jointly* phase small variants and SVs, pass both VCF files to HiPhase and specify two output VCF files in the same order as input.
Currently, DeepVariant and pbsv are the two supported input types.
While not required, it is *recommended* that global re-alignment is enabled when SVs are used.
If global re-alignment is used, a reference genome file _must_ be provided via `--reference {REFERENCE}`.

```bash
hiphase \
    --reference human_GRCh38_no_alt_analysis_set.fasta \
    --global-realignment-cputime 300 \
    --vcf HG001.GRCh38.deepvariant.vcf.gz \
    --output-vcf HG001.GRCh38.deepvariant.phased.vcf.gz \
    --vcf HG001.GRCh38.pbsv.vcf.gz \
    --output-vcf HG001.GRCh38.pbsv.phased.vcf.gz \
    ...
```

## Haplotagging
Haplotagging is the process of assigning a haplotype of origin to each read mapping, often for the purpose of visualizing data in IGV or for downstream analysis on haplotagged reads.
In HiPhase, haplotagging is done _while_ phasing (as opposed to a follow-up command).
To enable haplotagging, provide one output BAM file (`--output-bam`) for each input BAM file.
Note that enabling haplotagging will lead to longer run-times (though typically negligible increases in CPU time) due to increased file I/O from writing large BAM files.

```bash
hiphase \
    --bam m64109_200805_204709.GRCh38.bam \
    --output-bam m64109_200805_204709.GRCh38.haplotagged.bam \
    --bam m64109_200807_075817.GRCh38.bam \
    --output-bam m64109_200807_075817.GRCh38.haplotagged.bam \
    ...
```

## Input filtering
In some instances, filtering reads and/or variant calls by quality can improve the phasing quality, typically by reducing the number of switch errors generated in the outputs.
However, this can come at a cost of shortening phase block lengths due to the loss of heterozygous variants and/or reads that connect them.
The following options allow a user control over these parameters:
* `--min-mapq {MAPQ}` - Sets a minimum mapping quality (MAPQ) for reads to be used for solving the phasing problems. By default, this value is 5. Decreasing this value can lead to long run-times in noisier regions of the genome (e.g. centromeres), but may enable longer blocks as well. Note: reads with MAPQ below this threshold will still be haplotagged if that is enabled.
* `--min-vcf-qual {GQ}` - Sets a minimum genotype quality (GQ) for a heterozygous variant to be included in the phase solution. Note: any heterozygous variants with GQ less than this threshold will still be output as unphased variants.

## Multi-sample VCFs
By default, HiPhase will determine the sample name by looking at the first sample of the first VCF provided to it on the CLI.
If you have a multi-sample VCF file where the sample of interest is not the first one, you can override the default behavior by specifying the `--sample-name {SAMPLE}` parameter.
Note that in either case, the exact sample name _must_ be present in all provided VCF files or HiPhase will error out.
Additionally, HiPhase does not currently parse the sample names in the provided BAM files, so only provide BAM files containing reads for the sample of interest.

## Phasing singleton blocks
By default, HiPhase will effectively skip phase blocks that only contain a single variant, as these blocks are typically uninformative for downstream applications (i.e. the block will always generated a "0|1" solution for the single variant).
However, there may be situations where generating this information and/or haplotagging the reads is still useful.
To enable the phasing of these blocks, pass `--phase-singletons` to the CLI.
Note: this is likely to lead to increased compute / run time depending on the relative abundance of singletons in the dataset and the quantity of reads overlapping that data.

# Supported upstream processes
The following upstream processes are supported as inputs to HiPhase:

* Aligners (BAM files):
  * [pbmm2](https://github.com/PacificBiosciences/pbmm2) (recommended)
  * [minimap2](https://github.com/lh3/minimap2)
* Variant callers
  * [DeepVariant](https://github.com/google/deepvariant) - for SNV/indel
  * [pbsv](https://github.com/PacificBiosciences/pbsv) - for structural variants

Other upstream processes may work with HiPhase, but there is no official support for them at this time.

# Output files
## Phased VCF files
The primary output from HiPhase is one or more phased VCF files.
Following [VCF v4.2 specifications](https://samtools.github.io/hts-specs/VCFv4.2.pdf), each phased variant has the genotype field (GT) changed and an associated phase set tag (PS) added.
Example changes in the genotype field include: `0/1` -> `1|0` or `1/2` -> `1|2`.
Following recommendations of the VCF v4.2 specification, the PS tag is set to the position (`POS`, 1-based) of the first variant in the phase block.
Note that if multiple VCF files are provided as input (e.g. small variants and structural variants), the PS tags are shared.
This means that the first variant in the block may match a variant position from a _different_ VCF file (e.g. since there are fewer SVs, the PS tag is often referrring to a small variant).

HiPhase allows heterozygous variants to be converted to homozygous if that leads to a more optimal solution.
This is most commonly caused by noisy genomic regions and/or noisy indel calls.
When this happens, HiPhase _does not_ convert the output variant into the corresponding homozygous representation.
Instead, that variant is simply left as an unphased heterozygous call without a phase set ID (i.e. the VCF line is copied without changes).

## Haplotagged BAM files
HiPhase follows the same [haplotagging convention as WhatsHap](https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-haplotag), allowing for visualization tools (such as IGV) to use existing frameworks.
Each mapping is tagged with both a phase set ID (PS) and a haplotype ID (HP).
The phase set ID should match the values present in the corresponding phased VCF files.
The HP tag will be either a "1" or a "2" corresponding to the haplotypes in the phased VCF file.
Within a single phase block, all mappings with the same read name will have the same HP tag.
Mappings of the same read to _different_ phase blocks (e.g. to different chromosomes) are _not_ guaranteed to have matching HP tags.

In some instances, mappings within a block may not be used for solving the phasing problem (e.g. they are filtered).
These mappings will still be haplotagged as long as they overlap at least one phased variant in the final solution.
Mappings that do not overlap any phased variants or that equally / ambiguously match the two haplotypes will not be haplotagged.

## Auxiliary statistics files
### Phase block file (`--blocks-file`)
This CSV/TSV file contains information about the the phase blocks that were output by HiPhase.

Fields:
* `source_block_index` - The index of the phase problem within HiPhase. These values are 0-based and monotonically increasing. Duplicates can appear if an initial phase block was split while solving the phase problem (e.g. see block #10 in the example).
* `phase_block_id` - The phase block ID, which should match PS tags in both the VCF and BAM outputs
* `chrom` - the chromosome the block is on
* `start` - the position (POS) of the first variant in the block, 1-based
* `end` - the position (POS) of the last variant in the block, 1-based
* `num_variants` - the total number of phased variants in the phase block

Example:
```
source_block_index	phase_block_id	chrom	start	end	num_variants
0	10492	chr1	10492	31295	44
1	42043	chr1	42043	78515	7
2	97937	chr1	97937	204487	342
3	257717	chr1	257717	290912	85
4	352280	chr1	352280	382987	132
5	418014	chr1	418014	433277	4
6	495703	chr1	495703	523623	6
7	595856	chr1	595856	669651	131
8	688388	chr1	688388	814033	136
9	841166	chr1	841166	1150189	458
10	1169858	chr1	1169858	1362103	56
10	1382658	chr1	1382658	1424459	6
...
```

### Chromosome Summary File (`--summary-file`)
This CSV/TSV file contains chromosome-level summary statistics for all phase blocks on that chromosome.
Additionally, a chromosome labeled "all" contains aggregate statistics for all phase blocks generated by HiPhase.

Fields:
* `chromosome` - the chromosome label or "all", all other values are in reference to phase blocks only from this chromosome
* `num_variants` - the number of variants in the input VCF(s); only heterozygous or _non-reference_ homozygous calls are counted
* `num_heterozygous` - the number of heterozygous calls in the input VCF(s)
  * `num_heterozygous = num_phased + num_unphased`
* `num_phased` - the number of variant calls that were phased by HiPhase
* `num_unphased` - the number of heterozygous calls that were left unphased by HiPhase
* `num_het_snv` - the number of heterozygous single-nucleotide variants in the input VCF(s)
* `num_phased_snv` - the number of heterozygous single-nucleotide variants that were phased by HiPhase
* `num_blocks` - the number of distinct phase blocks generated, including `num_singletons`; this corresponds to the number of unique `PS` tags in the output
* `num_singletons` - the number of phase blocks with only one variant
* `variants_per_block_{function}` - the number of variants per block; the output is reduced by a function: median, mean, minimum (min), maximum (max), or summation (sum)
* `basepairs_per_block_{function}` - the number of basepairs per block; the output is reduced by a function: median, mean, minimum (min), maximum (max), or summation (sum)
* `block_ng50` - the largest phase block length such that all blocks of length greater than or equal to it span at least 50% of the chromosome

Example:
```
chromosome	num_variants	num_heterozygous	num_phased	num_unphased	num_het_snv	num_phased_snv	num_blocks	num_singletons	variants_per_block_median	variants_per_block_mean	variants_per_block_min	variants_per_block_max	variants_per_block_sum	basepairs_per_block_median	basepairs_per_block_mean	basepairs_per_block_min	basepairs_per_block_max	basepairs_per_block_sum	block_ng50
chr1	390387	227722	221689	6033	179630	176351	1118	169	11	198	1	3333	221689	45136	191472	1	18375626	214066410	526761
chr2	394773	239117	235123	3994	189802	188324	1181	190	8	199	1	4165	235123	40798	172873	1	2395966	204163229	449236
chr3	323127	198541	195171	3370	158502	157009	909	144	10	214	1	3369	195171	51392	186875	1	3576006	169869885	496188
chr4	349787	207079	203208	3871	166367	164539	910	153	8	223	1	4149	203208	41692	176245	1	3173866	160383760	463430
...
```

### Algorithm Statistics File (`--stats-file`)
This CSV/TSV file contains statistics regarding the performance of the underlying algorithms while running HiPhase.
This file is primarily for developers looking to improve HiPhase, but may be of use while identifying problematic phase blocks.
For details on how these relate to the internal algorithms, refer to [our methods](./methods.md).

Fields:
* `block_index` - the index of the block from phase block generation; this may correspond to multiple output phase blocks if the block was split; while these are monotonically increasing, they may show up out-of-order in the output file
* `chrom` - the chromosome for the block
* `start` - the position of the first variant in the block (0-based)
* `end` - the position of the last variant in the block (0-based)
* `num_variants` - the number of heterozygous variants in the block
* `num_reads` - the number of read mappings loaded to perform phasing
* `skipped_reads` - the number of read mappings that were loaded but not used for phasing
* `num_alleles` - the number of alleles that were defined (exact and in-exact) while loading reads
* `allele_matches` - the number of exactly-matching alleles for a given variant type; these are stored as an array corresponding to this order: [SNV, Insertion, Deletion, Indel, SvInsertion, SvDeletion, SvDuplication, SvInversion, SvBreakend, Unknown]
* `allele_partials` - the number of non-exact matching alleles for a given variant type (see order above); note that if global re-alignment is used, _all_ matching alleles will be assigned to this category
* `allele_failures` - the number of alleles that failed to match for a given variant type (see order above); this includes ambiguous / equal matches and deleted alleles
* `allele0_assigned` - the number of alleles internally assigned to allele0 (typically reference) for a given variant type (see order above)
* `allele1_assigned` - the number of alleles internally assigned to allele1 (always alternate) for a given variant type (see order above); for diploid organisms, `allele0_assigned` and `allele1_assigned` are expect to be roughly equal, and deviations may represent bias and/or lower quality variant inputs
* `is_global_realignment` - if `true`, then global re-alignment was used to extract alleles for this block and all stats reflect that process; otherwise, if `false`, then local re-alignment was used
* `pruned_solutions` - the number of solutions pruned during the A* traversal; if 0, then this block has a _guaranteed_ optimal solution given the problem definition
* `estimated_cost` - the estimated, heuristic cost to solve this block
* `actual_cost` - the actual cost of the solution from HiPhase
* `cost_ratio` - `estimated_cost / actual_cost`; ideally for a good heuristic, this value is near 1.0; blocks with significantly lower ratios indicates the heuristic did not accurately capture the cost of phasing the block, typically due to genomic complexities in the region
* `phased_variants` - the number of phased heterozygous variants in the HiPhase solution
* `homozygous_variants` - the number of heterozygous variants in the HiPhase solution that were converted to homozygous to achieve an optimal solution
* `skipped_variants` - the number of heterozygous variants that were skipped due to inability to phase it; this is usually caused by failure to match the allele on any overlapping mappings


Example:
```
block_index,chrom,start,end,num_variants,num_reads,skipped_reads,num_alleles,allele_matches,allele_partials,allele_failures,allele0_assigned,allele1_assigned,is_global_realignment,pruned_solutions,estimated_cost,actual_cost,cost_ratio,phased_variants,homozygous_variants,skipped_variants
1,chr1,42042,78514,7,15,21,60,"[10, 18, 0, 0, 0, 0, 0, 0, 0, 0]","[12, 18, 2, 0, 0, 0, 0, 0, 0, 0]","[0, 1, 0, 0, 0, 0, 0, 0, 0, 0]","[4, 6, 0, 0, 0, 0, 0, 0, 0, 0]","[18, 30, 2, 0, 0, 0, 0, 0, 0, 0]",false,0,11,11,1.0,7,0,0
0,chr1,10107,31294,62,76,183,2916,"[2239, 136, 22, 0, 0, 0, 0, 0, 0, 0]","[432, 85, 2, 0, 0, 0, 0, 0, 0, 0]","[13, 26, 0, 0, 0, 0, 0, 0, 0, 0]","[1492, 64, 19, 0, 0, 0, 0, 0, 0, 0]","[1179, 157, 5, 0, 0, 0, 0, 0, 0, 0]",false,3705,21294,23957,0.8888425094961807,44,18,0
...
```
