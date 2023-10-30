
/// A*-based phasing implementation
pub mod astar_phaser;
/// Functionality that iterates over VCF and BAM to form the prototype phase blocks
pub mod block_gen;
/// CLI functionality and checks
pub mod cli;
/// Contains multiple wrappers for useful data types in HiPhase
pub mod data_types;
/// Organizes primary workflow for a phase block including loading variants from VCF, loading reads from BAMs, running the phaser, and bundling the results
pub mod phaser;
/// Components for loading reads from a BAM file and converting them into haplotype observations
pub mod read_parsing;
/// Basic helpful utilities for pairwise sequence alignment
pub mod sequence_alignment;
/// Graph-based WFA - this is basically POA + WFA, but only allowing for measuring edit distance and no loops
pub mod wfa_graph;
/// Contains all the various output writer functionality
pub mod writers;
