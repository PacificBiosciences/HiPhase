# Methods
HiPhase is a tool for _jointly_ phasing multiple variant types with a single algorithm.
While there are algorithms that exist for phasing small variants and structural variants, most of them either 1) phase the variants sequentially (e.g. small variants first, and then overlay structural variants) or 2) phase only a subset of what HiPhase supports (e.g. small variants only).
To our knowledge, HiPhase is unique in that all supported variant types (see below) are phased simultaneously using a single algorithm.
As of v1.0.0, HiPhase supports the following variant types:

* Small variants from DeepVariant - SNVs, Insertions, Deletions, and Indels
* Structural variants (SVs) from pbsv - SvInsertions and SvDeletions
* Tandem repeats (TRs) from TRGT - TandemRepeat

In addition to the benefits of jointly phasing, there are also many advantages of using HiPhase over other phasing tools such as:

* [Long phase blocks with high accuracy](./performance.md)
* Gap spanning phase block generation - HiPhase includes logic to span coverage gaps caused by large deletions, inversions, and/or reference gaps
* Dual-mode allele assignment - for assigning alleles to reads, HiPhase implements two allele assignment methods
  * a local re-alignment algorithm - this is conceptually very similar to WhatsHap's allele assignment methods; this mode is recommended for phasing with only small variants
  * a global re-alignment algorithm - a novel application of the graph-based WFA algorithm to find an optimal allele assignment for a whole mapping; this mode is recommended for phasing with small variants, structural variants, and tandem repeats
* Novel phasing algorithm - HiPhase uses a novel application of the [A* search algorithm](https://en.wikipedia.org/wiki/A*_search_algorithm) to solve the phasing problem
* No downsampling - all provided read mappings are used; if the data has 100x coverage, HiPhase will use 100x coverage
* Multi-allelic variation is supported
* Quality of life additions:
  * Haplotagging while phasing - single command generates both phased VCF(s) and haplotagged BAM(s)
  * Built-in multi-threading - phase blocks are handled independently by different threads
  * Written in Rust - provides assurances on performance and memory safety, while also offering a framework for unit testing and distributing HiPhase

Greater detail on the internal HiPhase methods will be released in an upcoming pre-print.
