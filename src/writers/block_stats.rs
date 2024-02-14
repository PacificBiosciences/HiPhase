
use log::debug;
use rustc_hash::FxHashMap as HashMap;
use serde::Serialize;
use std::fs::File;
use std::path::Path;

use crate::block_gen::PhaseBlock;
use crate::data_types::reference_genome::ReferenceGenome;
use crate::data_types::variants::{VariantType, Zygosity};
use crate::phaser::PhaseResult;

/// This is a wrapper for writing out any stats to a file
#[derive(Default)]
pub struct BlockStatsCollector {
    /// Blocks that will be written out eventually
    blocks: Vec<PhaseBlock>,
    /// Tracks the number of phased SNVs, key is (sample_name, chromosome)
    phased_snvs: HashMap<(String, String), usize>
}

/// Contains all the data written to each row of our stats file
#[derive(Serialize)]
struct BlockRow {
    /// The index of the block
    source_block_index: usize,
    /// the sample name tied to the block
    sample_name: String,
    /// Phase set ID - usually the position first variant in the block
    phase_block_id: u64,
    /// the chromosome of the block
    chrom: String,
    /// the position of the first allele
    start: u64,
    /// the position of the last allele
    end: u64,
    /// the number of variants in the block
    num_variants: usize
}

#[derive(Serialize)]
struct SummaryRow {
    /// the sample name
    sample_name: String,
    /// The chromosome or "all"
    chromosome: String,
    /// The total number of variants
    num_variants: usize,
    /// The total number of heterozygous variants
    num_heterozygous: usize,
    /// The number of phased heterozygous variants
    num_phased: usize,
    /// The number of unphased heterozygous variants
    num_unphased: usize,
    /// The number of heterozygous SNVs
    num_het_snv: usize,
    /// The number of phased, heterozygous SNVs
    num_phased_snv: usize,
    /// The total number of blocks
    num_blocks: usize,
    /// The number of blocks with only 1 phased variant
    num_singletons: usize,
    /// variants per block stats
    variants_per_block_median: usize,
    variants_per_block_mean: usize,
    variants_per_block_min: usize,
    variants_per_block_max: usize,
    variants_per_block_sum: usize,
    /// basepairs per block stats
    basepairs_per_block_median: u64,
    basepairs_per_block_mean: u64,
    basepairs_per_block_min: u64,
    basepairs_per_block_max: u64,
    basepairs_per_block_sum: u64,
    /// block NG50
    block_ng50: Option<u64>
}

impl BlockStatsCollector {
    /// Creates a new writer for a given filename
    /// # Arguments
    /// * `csv_filename` - the path to write all stats to
    pub fn new() -> BlockStatsCollector {
        Self::default()
    }

    /// Adds a block to our collection
    /// # Arguments
    /// * `block` - the block to add, no checks are performed on the input
    pub fn add_block(&mut self, block: PhaseBlock) {
        self.blocks.push(block);
    }

    /// Adds a phase result to our statistics
    /// # Arguments
    /// * `chrom` - the chromosome for the result
    /// * `result` 
    pub fn add_result(&mut self, result: &PhaseResult) {
        if let Some(stats) = result.statistics.as_ref() {
            if let Some(count) = stats.phased_snvs() {
                let sample_name: String = result.phase_block.sample_name().to_string();
                let chrom: String = result.phase_block.get_chrom().to_string();
                *self.phased_snvs.entry((sample_name, chrom)).or_insert(0) += count as usize;
            }
        }
    }

    /// Will write all blocks to a CSV filename in order
    /// # Arguments
    /// * `filename` - the filename for the output (tsv/csv)
    pub fn write_blocks(&mut self, filename: &Path) -> csv::Result<()> {
        // modify the delimiter to "," if it ends with .csv
        let is_csv: bool = filename.extension().unwrap_or_default() == "csv";
        let delimiter: u8 = if is_csv { b',' } else { b'\t' };
        let mut csv_writer: csv::Writer<File> = csv::WriterBuilder::new()
            .delimiter(delimiter)
            .from_path(filename)?;

        // make sure we go through the blocks in order
        self.blocks.sort();
        for block in self.blocks.iter() {
            let block_row = BlockRow {
                source_block_index: block.get_block_index(),
                sample_name: block.sample_name().to_string(),
                phase_block_id: block.get_start()+1,
                chrom: block.get_chrom().to_string(),
                start: block.get_start()+1,
                end: block.get_end()+1,
                num_variants: block.get_num_variants()
            };
            csv_writer.serialize(&block_row)?;
        }
        csv_writer.flush()?;
        Ok(())
    }

    /// Will write out a file containing chromosome level block statistics and overall statistics.
    /// # Arguments
    /// * `filename` - the filename for the output (tsv/csv)
    /// * `reference_genome` - the reference genome, without this we can't determine lengths
    /// * `variant_counts` - the variant count data from parsing; key is (sample, chromosome, variant type, zygosity), value is a count
    pub fn write_block_stats(
        &self, 
        sample_order: &[String],
        filename: &Path, 
        reference_genome: &ReferenceGenome, 
        variant_counts: HashMap<(String, String, VariantType, Zygosity), usize>
    ) -> csv::Result<()> {
        // modify the delimiter to "," if it ends with .csv
        let is_csv: bool = filename.extension().unwrap_or_default() == "csv";
        let delimiter: u8 = if is_csv { b',' } else { b'\t' };
        let mut csv_writer: csv::Writer<File> = csv::WriterBuilder::new()
            .delimiter(delimiter)
            .from_path(filename)?;
        
        // calculate the total reference length assuming we have a reference
        let total_contig_length: u64 = {
            let mut contig_sum: u64 = 0;
            for contig_name in reference_genome.contig_keys().iter() {
                let contig_length = reference_genome.get_full_chromosome(contig_name).len();
                contig_sum += contig_length as u64;
            }
            contig_sum
        };

        for sample_name in sample_order.iter() {
            // go through all the blocks and generate chromosome level and overall stats
            let mut blocks_by_chrom: HashMap<String, Vec<PhaseBlock>> = Default::default();
            let mut all_sample_blocks: Vec<PhaseBlock> = Default::default();
            
            for block in self.blocks.iter() {
                if block.sample_name() == sample_name {
                    let chrom: String = block.get_chrom().to_string();
                    blocks_by_chrom.entry(chrom).or_default().push(block.clone());
                    all_sample_blocks.push(block.clone());
                }
            }

            // these are from parsing the VCFs
            let mut num_variants: HashMap<String, usize> = Default::default();
            let mut num_heterozygous: HashMap<String, usize> = Default::default();
            let mut num_het_snv: HashMap<String, usize> = Default::default();
            let mut ordered_iteration: Vec<_> = variant_counts.iter().collect();
            ordered_iteration.sort();
            for ((sample, chrom, variant_type, zygosity), &count) in ordered_iteration.iter() {
                debug!("{} {} {:?} {:?} => {}", sample, chrom, variant_type, zygosity, count);
                if sample == sample_name && *variant_type != VariantType::Unknown && *zygosity != Zygosity::HomozygousReference && *zygosity != Zygosity::Unknown {
                    // okay, we're working with something we can count now
                    *num_variants.entry(chrom.clone()).or_insert(0) += count;

                    if *zygosity == Zygosity::Heterozygous {
                        *num_heterozygous.entry(chrom.clone()).or_insert(0) += count;
                        if *variant_type == VariantType::Snv {
                            *num_het_snv.entry(chrom.clone()).or_insert(0) += count;
                        }
                    }
                }
            }

            // best place to pull contig order from is our reference genome
            for contig in reference_genome.contig_keys().iter() {
                let contig_length: u64 = reference_genome.get_full_chromosome(contig).len() as u64;
                let chrow_stats_row: SummaryRow = self.generate_summary_row(
                    sample_name,
                    contig, 
                    blocks_by_chrom.get(contig).unwrap_or(&vec![]), 
                    *num_variants.get(contig).unwrap_or(&0), 
                    *num_heterozygous.get(contig).unwrap_or(&0), 
                    *num_het_snv.get(contig).unwrap_or(&0),
                    *self.phased_snvs.get(&(sample_name.clone(), contig.clone())).unwrap_or(&0),
                    contig_length
                );
                csv_writer.serialize(&chrow_stats_row)?;
            }

            let all_stats_row: SummaryRow = self.generate_summary_row(
                sample_name, 
                "all",
                &all_sample_blocks,
                num_variants.values().sum(),
                num_heterozygous.values().sum(),
                num_het_snv.values().sum(),
                self.phased_snvs.iter().map(|(_chrom, &count)| count).sum(),
                total_contig_length
            );
            csv_writer.serialize(&all_stats_row)?;
        }

        csv_writer.flush()?;
        Ok(())
    }

    /// Utility function for building a summary row from given data.
    /// # Arguments
    /// * `sample_name` - pass-through sample ID
    /// * `chrom` - pass-through chromosome name
    /// * `blocks` - the phase blocks that get parsed into a summary row for the sample/chromosome
    /// * `num_variants` - total number of variants
    /// * `num_heterozygous` - total number of heterozygous variants
    /// * `num_het_snv` - total number of heterozygous SNVs
    /// * `num_phased_snv` - total number of output phased SNVs
    /// * `contig_length` - length of the chromosome
    #[allow(clippy::too_many_arguments)]
    fn generate_summary_row(
        &self, sample_name: &str, chrom: &str, blocks: &[PhaseBlock],
        num_variants: usize, num_heterozygous: usize, num_het_snv: usize, num_phased_snv: usize, 
        contig_length: u64
    ) -> SummaryRow {
        // make sure every block has the correct sample name
        assert!(blocks.iter().all(|b| b.sample_name() == sample_name));

        // these are derivable from our results
        let num_blocks = blocks.len();
        let num_singletons = blocks.iter()
            .filter(|block| block.get_num_variants() == 1)
            .count();

        // collect the block variant stats, this is only counting heterozygous variants
        let mut block_variants: Vec<usize> = blocks.iter()
            .map(|b| b.get_num_variants())
            .collect();
        
        // collect the block length stats
        let mut block_lengths: Vec<u64> = blocks.iter()
            .map(|b| b.bp_len())
            .collect();

        // let total_heterozygous: usize = num_heterozygous.iter().map(|(_chrom, &count)| count).sum();
        let num_phased: usize = block_variants.iter().sum();
        let num_unphased = num_heterozygous - num_phased; 

        block_variants.sort();
        let variants_per_block_median: usize = if block_variants.is_empty() { 0 } else { block_variants[block_variants.len() / 2] };
        let variants_per_block_mean: usize = if block_variants.is_empty() { 0 } else { block_variants.iter().sum::<usize>() / block_variants.len() };
        let variants_per_block_min: usize = *block_variants.iter().min().unwrap_or(&0);
        let variants_per_block_max: usize = *block_variants.iter().max().unwrap_or(&0);
        let variants_per_block_sum: usize = block_variants.iter().sum();

        block_lengths.sort();
        let basepairs_per_block_median: u64 = if block_lengths.is_empty() { 0 } else { block_lengths[block_lengths.len() / 2] };
        let basepairs_per_block_mean: u64 = if block_lengths.is_empty() { 0 } else { block_lengths.iter().sum::<u64>() / block_lengths.len() as u64 };
        let basepairs_per_block_min: u64 = *block_lengths.iter().min().unwrap_or(&0);
        let basepairs_per_block_max: u64 = *block_lengths.iter().max().unwrap_or(&0);
        let basepairs_per_block_sum: u64 = block_lengths.iter().sum();

        let block_ng50: Option<u64> = if contig_length != 0 {
            Some(calculate_block_ng50(&block_lengths, contig_length))
        } else {
            None
        };

        SummaryRow {
            sample_name: sample_name.to_string(),
            chromosome: chrom.to_string(),
            num_variants,
            num_heterozygous,
            num_phased,
            num_unphased,
            num_het_snv,
            num_phased_snv,
            num_blocks,
            num_singletons,
            variants_per_block_median,
            variants_per_block_mean,
            variants_per_block_min,
            variants_per_block_max,
            variants_per_block_sum,
            basepairs_per_block_median,
            basepairs_per_block_mean,
            basepairs_per_block_min,
            basepairs_per_block_max,
            basepairs_per_block_sum,
            block_ng50
        }
    }
}

/// Helper subroutine for calculating block NG50 from a list of sorted blocks.
/// # Arguments
/// * `sorted_blocks` - block sizes sorted in ascending order
/// * `contig_length` - the total contig length, half is needed to reach NG50
/// # Panics
/// * if while iterating it detects unsorted blocks
fn calculate_block_ng50(sorted_blocks: &[u64], contig_length: u64) -> u64 {
    let mut last_block_size: u64 = u64::MAX;
    let mut length_sum: u64 = 0;

    // add one to handle odd values (e.g. rounding up)
    let target_length: u64 = (contig_length + 1) / 2;

    for &block_size in sorted_blocks.iter().rev() {
        // we're going in reverse, so block sizes *should* be monotonically decreasing
        assert!(block_size <= last_block_size);
        last_block_size = block_size;

        // add in the block and check again half the total length
        length_sum += block_size;
        if length_sum >= target_length {
            // we made it, return the size of this block
            return block_size;
        }
    }

    // we didn't make it
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculate_block_ng50() {
        // odd block length
        let contig_length: u64 = 21;
        let blocks: Vec<u64> = vec![1, 2, 3, 4, 10];
        let bad_blocks: Vec<u64> = vec![2];
        let good_blocks: Vec<u64> = vec![9, 10];

        assert_eq!(calculate_block_ng50(&blocks, contig_length), 4);
        assert_eq!(calculate_block_ng50(&bad_blocks, contig_length), 0);
        assert_eq!(calculate_block_ng50(&good_blocks, contig_length), 9);

        // same but now even block length
        let contig_length: u64 = 20;
        assert_eq!(calculate_block_ng50(&blocks, contig_length), 10);
        assert_eq!(calculate_block_ng50(&bad_blocks, contig_length), 0);
        assert_eq!(calculate_block_ng50(&good_blocks, contig_length), 10);
    }
}