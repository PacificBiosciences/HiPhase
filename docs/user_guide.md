# User Guide
Table of contents:

* [Quickstart](#quickstart)
* [Common uses cases](#common-use-cases)
* [Supported upstream processes](#supported-upstream-processes)
* [Output files](#output-files)
* [FAQ](#faq)

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
* `--bam {IN_BAM}` - path to a BAM file containing reads only from the sample that is being phased, this option can be specified multiple times; each BAM file can only contain reads from a single sample; each BAM file must be indexed prior to running HiPhase
* `--vcf {IN_VCF}` - path to a VCF file containing the variants to phase, this option can be specified multiple times; each sample being phased must appear in each provided VCF file; we recommend providing at least a small variant VCF and structural variant VCF; each VCF file must be index prior to running HiPhase
* `--output-vcf {OUT_VCF}` - path to the output VCF that will contain the phased variants, this option must be specified the same number of times as `--vcf` 
* `--reference {REFERENCE}` - a FASTA file containing the reference genome, gzip allowed; as of v0.9.0, this is a required parameter
* `--threads {THREADS}` - number of threads to use for phasing (default: 1)

## Quickstart Example
This example also includes some auxiliary outputs that may be useful for tracking phase result statistics.
This example also only uses a small variant VCF, so `--disable-global-realignment` is recommended.
For details on each, refer to the [Output files](#output-files) section.

```
hiphase \
    --disable-global-realignment \
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

## Recommended resources
HiPhase has built in parallel processing via the `--threads` parameter.
We recommend reserving 4 GB of memory per thread allocated to HiPhase.
For example, most of our internal tests use 16 threads and reserve 64 GB of memory.

# Common use cases
## Joint phasing small variants, structural variants, and tandem repeats
To *jointly* phase small variants, structural variants, and tandem repeats, pass all VCF files to HiPhase and specify one output VCF file for each in the same order as input.
Currently, DeepVariant, pbsv, and TRGT are the three supported input types.
Global re-alignment is *strongly recommended* (and enabled by default) when structural variants and tandem repeats are provided.

```bash
hiphase \
    --reference human_GRCh38_no_alt_analysis_set.fasta \
    --vcf HG001.GRCh38.deepvariant.vcf.gz \
    --output-vcf HG001.GRCh38.deepvariant.phased.vcf.gz \
    --vcf HG001.GRCh38.pbsv.vcf.gz \
    --output-vcf HG001.GRCh38.pbsv.phased.vcf.gz \
    --vcf HG001.GRCh38.trgt.vcf.gz \
    --output-vcf HG001.GRCh38.trgt.phased.vcf.gz \
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

To phase multiple samples at once from a single VCF file, each sample must be specified via the `--sample-name` parameter.
Additionally, all BAM files for the samples will need to be provided with the `--bam` parameter.
These files can be haplotagged as normal via the `--output-bam` option.
Currently, multi-sample BAM files are _not_ supported, so each provided BAM *must* match with exactly one VCF sample identifier according to the read group (RG) sample name (SM) tag.

## Phasing singleton blocks
By default, HiPhase will effectively skip phase blocks that only contain a single variant, as these blocks are typically uninformative for downstream applications (i.e. the block will always generated a "0|1" solution for the single variant).
However, there may be situations where generating this information and/or haplotagging the reads is still useful.
To enable the phasing of these blocks, pass `--phase-singletons` to the CLI.
Note: this is likely to lead to increased compute / run time depending on the relative abundance of singletons in the dataset and the quantity of reads overlapping that data.

## Controlling dual-mode allele assignment
In HiPhase, dual-mode allele assignment is enabled by default.
This process attempts to globally realign a mapping to assign alleles, reverting to local realignment as a backup.
In general, global realignment is more accurate, but also more computationally costly in noisy genomic regions and/or regions with incomplete variant calls.
To prevent excessive run-times, there are checks in place that will stop global realignment for an individual mapping and revert to local realignment if it becomes too expensive.
Additionally, if too many mappings within a block are failing global realignment (>50% by default), then HiPhase will revert the rest of the block to using local realignment only.
This combination of checks allows HiPhase to switch between global and local realignment in a deterministic way.
The following options control how global realignment is used for read parsing:

* `--disable-global-realignment` - this will disable global realignment for the entire process; this is recommended if _only_ a small variant VCF is available
* `--global-realignment-max-ed <DISTANCE>` - sets a maximum on the edit distance for global realignment for a single read mapping; if this maximum is reached, HiPhase will revert to local realignment for the mapping; increasing this value may increase run-time but also generate better realignments
* `--max-global-failure-ratio <FRAC>` - sets a maximum fraction of global realignment failures to allow before reverting the rest of the phase block to local realignment; increasing this may lead to more global realignments, but also to longer run times; decreasing this may lead to fewer global realignments, but also improved run times in noisy phase blocks
* `--global-failure-count <COUNT>` - sets the minimum number of global realignment failures before the failure ratio check is active 

# Supported upstream processes
The following upstream processes are supported as inputs to HiPhase:

* Aligners (BAM files):
  * [pbmm2](https://github.com/PacificBiosciences/pbmm2) (recommended)
  * [minimap2](https://github.com/lh3/minimap2)
* Variant callers
  * [DeepVariant](https://github.com/google/deepvariant) - for SNV/indel
  * [pbsv](https://github.com/PacificBiosciences/pbsv) - for structural variants
  * [TRGT](https://github.com/PacificBiosciences/trgt) - for tandem repeats

Other upstream processes may work with HiPhase, but there is no official support for them at this time.

# Output files
## Phased VCF files
The primary output from HiPhase is one or more phased VCF files.
Following [VCF v4.2 specifications](https://samtools.github.io/hts-specs/VCFv4.2.pdf), each phased variant has the genotype field (GT) changed and an associated phase set tag (PS) added.
Example changes in the genotype field include: `0/1` -> `1|0` or `1/2` -> `1|2`.
Following recommendations of the VCF v4.2 specification, the PS tag is set to the position (`POS`, 1-based) of the first variant in the phase block.
Note that if multiple VCF files are provided as input (e.g. small variants and structural variants), the PS tags are shared.
This means that the first variant in the block may match a variant position from a _different_ VCF file (e.g. since there are fewer SVs, the PS tag is often referrring to a small variant).

HiPhase allows heterozygous variants to be treated as homozygous if that leads to a more optimal solution.
This is most commonly caused by noisy genomic regions and/or noisy indel calls.
When this happens, HiPhase _does not_ convert the output variant into the corresponding homozygous representation.
Instead, that variant is simply left as an unphased heterozygous call without a phase set ID (i.e. the VCF line is copied without changes).

### Additional VCF fields
The following are non-standard tags that HiPhase may include in the VCF output:

* FORMAT field `PF` tag - Short for **P**hase **F**lags. This will be set with a String label when the variant was intentionally excluded from the phasing algorithm. Possible values:
  * `TR_OVERLAP` - indicates that the variant is fully contained within a tandem repeat call and was intentionally ignored because it is likely less complete and/or correct when compared to the tandem repeat call

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
* `source_block_index` - The index of the phase problem within HiPhase. These values are 0-based and monotonically increasing. Duplicates can appear if an initial phase block was split while solving the phase problem (e.g. see block #10 in the example). Blocks indices may be skipped if singleton blocks are not enabled.
* `sample_name` - the sample name for the block, mostly for multi-sample inputs
* `phase_block_id` - The phase block ID, which should match PS tags in both the VCF and BAM outputs
* `chrom` - the chromosome the block is on
* `start` - the position (POS) of the first variant in the block, 1-based
* `end` - the position (POS) of the last variant in the block, 1-based
* `num_variants` - the total number of phased variants in the phase block

Example:
```
source_block_index	sample_name	phase_block_id	chrom	start	end	num_variants
0	HG001	10492	chr1	10492	31295	58
1	HG001	42043	chr1	42043	78515	7
2	HG001	97937	chr1	97937	157660	131
3	HG001	166982	chr1	166982	204487	211
4	HG001	257717	chr1	257717	290912	91
5	HG001	352280	chr1	352280	382987	153
6	HG001	418014	chr1	418014	433277	4
7	HG001	495703	chr1	495703	523623	6
8	HG001	595856	chr1	595856	669651	131
...
```

### Chromosome Summary File (`--summary-file`)
This CSV/TSV file contains chromosome-level summary statistics for all phase blocks on that chromosome.
Additionally, a chromosome labeled "all" contains aggregate statistics for all phase blocks generated by HiPhase.

Fields:
* `sample_name` - the sample name for the statistics, mostly for multi-sample VCF inputs
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
sample_name	chromosome	num_variants	num_heterozygous	num_phased	num_unphased	num_het_snv	num_phased_snv	num_blocks	num_singletons	variants_per_block_median	variants_per_block_mean	variants_per_block_min	variants_per_block_max	variants_per_block_sum	basepairs_per_block_median	basepairs_per_block_mean	basepairs_per_block_min	basepairs_per_block_max	basepairs_per_block_sum	block_ng50
HG001	chr1	390387	227722	222547	5175	179630	177207	1022	55	15	217	1	3333	222547	57130	207513	1	18375626	212078774	481758
HG001	chr2	394773	239117	235274	3843	189802	188509	1063	63	12	221	1	4165	235274	54325	190369	1	2395966	202362848	437966
HG001	chr3	323127	198541	195312	3229	158502	157181	818	44	16	238	1	3369	195312	68954	205458	1	3374952	168064953	482656
HG001	chr4	349787	207079	203594	3485	166367	164916	823	61	12	247	1	4148	203594	52601	194864	1	3173866	160373348	463430
...
```

### Haplotag file (`--haplotag-file`)
This CSV/TSV file contains haplotag information for aligned reads. 
Note that while this contains the same information as the HP tag in the haplotagged output BAMs, generating those output BAMs is not required to generate this file.

Fields:
* `source_block_index` - The index of the phase problem within HiPhase. These values are 0-based and monotonically increasing. Duplicates can appear if an initial phase block was split while solving the phase problem. Blocks indices may be skipped if singleton blocks are not enabled.
* `sample_name` - The sample name for the block, mostly for multi-sample inputs.
* `chrom` - The chromosome the block is on.
* `phase_block_id` - The phase block ID, which should match PS tags in both the VCF and BAM outputs.
* `read_name` - The read name from the BAM file, all mappings for a given read will have the same haplotag within a single block.
* `haplotag` - The assigned haplotype ID (HP in BAM), either 1 or 2.

Example:
```
source_block_index	sample_name	chrom	phase_block_id	read_name	haplotag
0	HG001	chr1	10492	m64109_200807_075817/118227363/ccs	1
0	HG001	chr1	10492	m64109_200810_062248/32113861/ccs	2
0	HG001	chr1	10492	m64109_200807_075817/6948622/ccs	1
0	HG001	chr1	10492	m64109_200813_162416/12846754/ccs	1
0	HG001	chr1	10492	m64109_200805_204709/170461776/ccs	1
0	HG001	chr1	10492	m64109_200815_033514/155779891/ccs	2
...
```

### Algorithm Statistics File (`--stats-file`)
This CSV/TSV file contains statistics regarding the performance of the underlying algorithms while running HiPhase.
This file is primarily for developers looking to improve HiPhase, but may be of use while identifying problematic phase blocks.
For details on how these relate to the internal algorithms, refer to [our methods](./methods.md).

Fields:
* `block_index` - the index of the block from phase block generation; this may correspond to multiple output phase blocks if the block was split; while these are monotonically increasing, they may show up out-of-order in the output file
* `sample_name` - the sample name for the block, mostly for multi-sample VCF inputs
* `chrom` - the chromosome for the block
* `start` - the position of the first variant in the block (0-based)
* `end` - the position of the last variant in the block (0-based)
* `num_variants` - the number of heterozygous variants in the block
* `num_reads` - the number of read mappings loaded to perform phasing
* `skipped_reads` - the number of read mappings that were loaded but not used for phasing
* `num_alleles` - the number of alleles that were defined (exact and in-exact) while loading reads
* `allele_matches` - the number of exactly-matching alleles for a given variant type; these are stored as an array corresponding to this order: [SNV, Insertion, Deletion, Indel, SvInsertion, SvDeletion, SvDuplication, SvInversion, SvBreakend, TandemRepeat, Unknown]
* `allele_partials` - the number of non-exact matching alleles for a given variant type (see order above); note that if global re-alignment is used, _all_ matching alleles will be assigned to this category
* `allele_failures` - the number of alleles that failed to match for a given variant type (see order above); this includes ambiguous / equal matches and deleted alleles
* `allele0_assigned` - the number of alleles internally assigned to allele0 (typically reference) for a given variant type (see order above)
* `allele1_assigned` - the number of alleles internally assigned to allele1 (always alternate) for a given variant type (see order above); for diploid organisms, `allele0_assigned` and `allele1_assigned` are expect to be roughly equal, and deviations may represent bias and/or lower quality variant inputs
* `global_aligned` - the number of mappings that were globally realigned
* `local_aligned` - the number of mappings that were locally realigned
* `pruned_solutions` - the number of solutions pruned during the A* traversal; if 0, then this block has a _guaranteed_ optimal solution given the problem definition
* `estimated_cost` - the estimated, heuristic cost to solve this block
* `actual_cost` - the actual cost of the solution from HiPhase
* `cost_ratio` - `estimated_cost / actual_cost`; ideally for a good heuristic, this value is near 1.0; blocks with significantly lower ratios indicates the heuristic did not accurately capture the cost of phasing the block, typically due to genomic complexities in the region
* `phased_variants` - the number of phased heterozygous variants in the HiPhase solution
* `homozygous_variants` - the number of heterozygous variants in the HiPhase solution that were converted to homozygous to achieve an optimal solution
* `skipped_variants` - the number of heterozygous variants that were skipped due to inability to phase it; this is usually caused by failure to match the allele on any overlapping mappings


Example:
```
block_index,sample_name,chrom,start,end,num_variants,num_reads,skipped_reads,num_alleles,allele_matches,allele_partials,allele_failures,allele0_assigned,allele1_assigned,global_aligned,local_aligned,pruned_solutions,estimated_cost,actual_cost,cost_ratio,phased_variants,homozygous_variants,skipped_variants
0,HG001,chr1,10107,31294,63,102,69,1779,"[10, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]","[1500, 148, 11, 0, 0, 109, 0, 0, 0, 0, 0]","[7, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0]","[878, 38, 6, 0, 0, 103, 0, 0, 0, 0, 0]","[632, 111, 5, 0, 0, 6, 0, 0, 0, 0, 0]",149,2,0,18040,18080,0.9977876106194691,58,5,0
1,HG001,chr1,42042,78514,8,10,9,40,"[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]","[14, 19, 2, 0, 0, 0, 0, 0, 0, 5, 0]","[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]","[12, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0]","[2, 19, 2, 0, 0, 0, 0, 0, 0, 0, 0]",18,0,0,0,0,1.0,8,0,0
...
```

# FAQ
## How do I fix "Error during BAM read group parsing: BAM file has no read groups (RG) tag"?
By default, HiPhase checks read group IDs to assign BAM files to a VCF sample ID for phasing.
If you are sure that all provided BAM files all correspond to a single sample ID, you can pass the `--ignore-read-groups` flag to disable this check.
Note that this mode can only be used for single-sample phasing.

## Why are some small/structural variants unphased when I added tandem repeats?
DeepVariant, pbsv, and TRGT all have overlap in the variants that are called, but this is especially pronounced with TRGT calls.
In general, the tandem repeat calls from TRGT are more accurate and less fragmented than the corresponding small variants or structural variant calls, which tend to have more errors or spurious calls and tend to get fragmented into multiple entries in the VCF file.
HiPhase handles these overlaps by ignoring (i.e., not phasing) any variant calls that are fully contained within a TRGT tandem repeat call.
Since TRGT calls often correspond to multiple small variant calls, the total number of phased variants tends to drop when TRGT is added.
Additionally, this has the benefit of removing false heterozygous variants in the tandem repeat regions.
This can have the side effect of (correctly) reducing block NG50 when those false variants are the only links between two otherwise unlinked phase blocks.
Variants that are intentionally excluded from phasing this way will have a FORMAT `PF` tag of `TR_OVERLAP`.
