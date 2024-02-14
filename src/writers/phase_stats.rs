
use crate::data_types::variants::VariantType;
use crate::phaser::PhaseResult;

use serde::Serialize;
use std::fs::File;
use std::ops::AddAssign;
use std::path::Path;

/// Contains statistics on the loading and parsing of reads into alleles.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadStats {
    /// The number of reads loaded
    num_reads: u64,
    /// The number of reads without any determined alleles
    skipped_reads: u64,
    /// The number of alleles successfully loaded
    num_alleles: u64,
    /// Records the number of exact matches found for a type
    exact_matches: [u64; VariantType::Unknown as usize + 1],
    /// Records the number of inexact matches found for a type
    inexact_matches: [u64; VariantType::Unknown as usize + 1],
    /// Records the number of no-matches found for a type
    failed_matches: [u64; VariantType::Unknown as usize + 1],
    /// Records the number of matches to allele 0
    allele0_matches: [u64; VariantType::Unknown as usize + 1],
    /// Records the number of matches to allele 1
    allele1_matches: [u64; VariantType::Unknown as usize + 1],
    /// The number of reads that used global re-alignment
    global_aligned: usize,
    /// The number of reads that used local re-alignment
    local_aligned: usize
}

impl ReadStats {
    /// Creates a new `ReadStats` struct and does some sanity checks.
    /// # Arguments
    /// * `num_reads` - the number of reads loaded
    /// * `skipped_reads` - the number of reads that were skipped during loading
    /// * `num_alleles` - the number of alleles successfully (i.e. un-ambiguously) loaded
    /// * `exact_matches` - the number of exact matches for each type
    /// * `inexact_matches` - the number of inexact matches for each type
    /// * `failed_matches` - the number of ambiguous or no-sequence alleles
    /// * `allele0_matches` - the number of alleles assigned to allele 0
    /// * `allele1_matches` - the number of alleles assigned to allele 1
    /// * `global_aligned` - the number of reads that were globally aligned
    /// * `local_aligned` - the number of reads that were locally aligned
    /// # Panics
    /// * if `num_reads` > `num_alleles`, because that would imply some reads have no alleles
    /// * if `num_alleles != exact_matches.sum() + inexact_matches.sum()`, because that would imply some alleles are not being counted correctly somewhere
    /// * if `num_alleles != allele0_matches.sum() + allele1_matches.sum()`, because that would imply some alleles are not being counted correctly somewhere
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        num_reads: u64, skipped_reads: u64, num_alleles: u64, 
        exact_matches: [u64; VariantType::Unknown as usize + 1], 
        inexact_matches: [u64; VariantType::Unknown as usize + 1], 
        failed_matches: [u64; VariantType::Unknown as usize + 1], 
        allele0_matches: [u64; VariantType::Unknown as usize + 1], 
        allele1_matches: [u64; VariantType::Unknown as usize + 1],
        global_aligned: usize,
        local_aligned: usize
    ) -> ReadStats {
        assert!(num_alleles >= num_reads);
        assert_eq!(num_alleles, exact_matches.iter().sum::<u64>() + inexact_matches.iter().sum::<u64>());
        assert_eq!(num_alleles, allele0_matches.iter().sum::<u64>() + allele1_matches.iter().sum::<u64>());
        ReadStats {
            num_reads,
            skipped_reads,
            num_alleles,
            exact_matches,
            inexact_matches,
            failed_matches,
            allele0_matches,
            allele1_matches,
            global_aligned,
            local_aligned
        }
    }
    
    // adders
    pub fn increase_num_reads(&mut self, value: u64) {
        self.num_reads += value;
    }

    pub fn increase_skipped_reads(&mut self, value: u64) {
        self.skipped_reads += value;
    }

    // getters
    pub fn skipped_reads(&self) -> u64 {
        self.skipped_reads
    }

    pub fn local_aligned(&self) -> usize {
        self.local_aligned
    }

    pub fn global_aligned(&self) -> usize {
        self.global_aligned
    }

    pub fn total_aligned(&self) -> usize {
        self.local_aligned + self.global_aligned
    }
}

impl AddAssign for ReadStats {
    fn add_assign(&mut self, rhs: Self) {
        self.num_reads += rhs.num_reads;
        self.skipped_reads += rhs.skipped_reads;
        self.num_alleles += rhs.num_alleles;
        add_assign_array(&mut self.exact_matches, &rhs.exact_matches);
        add_assign_array(&mut self.inexact_matches, &rhs.inexact_matches);
        add_assign_array(&mut self.failed_matches, &rhs.failed_matches);
        add_assign_array(&mut self.allele0_matches, &rhs.allele0_matches);
        add_assign_array(&mut self.allele1_matches, &rhs.allele1_matches);
        self.global_aligned += rhs.global_aligned;
        self.local_aligned += rhs.local_aligned;
    }
}

/// Wrapper function for adding two array of equal length element wise
fn add_assign_array(lhs: &mut [u64], rhs: &[u64]) {
    assert_eq!(lhs.len(), rhs.len());
    for (l, &r) in lhs.iter_mut().zip(rhs.iter()) {
        *l += r;
    }
}

/// Contains any statistics from the phasing problem solver that may be relevant
pub struct PhaseStats {
    /// The number of solutions that were pruned during calculation
    pruned_solutions: Option<u64>,
    /// For heuristic solvers, the estimate cost prior to computation
    estimated_cost: Option<u64>,
    /// The actual cost of the solution
    actual_cost: Option<u64>,
    /// The number of phased variants
    phased_variants: Option<u64>,
    /// The number of phased SNV variants
    phased_snvs: Option<u64>,
    /// The number of variants where the phasing solution turned them homozygous
    homozygous_variants: Option<u64>,
    /// The number of ignored variants
    skipped_variants: Option<u64>
}

impl PhaseStats {
    /// Creates phase stats for an A* algorithm solution
    /// # Arguments
    /// * `pruned_solutions` - the number of solutions that were pruned from the solution space to reduce memory consumption and run-time; if this is 0, then we have a guaranteed best solution
    /// * `estimated_cost` - the estimated cost of this phase block calculated from the heuristic
    /// * `actual_cost` - the actual cost of the final solution
    /// * `phased_variants` - the number of variants that were phased in the solution (e.g. not converted to homozygous)
    /// * `homozygous_variants` - the number of variants that were converted to homozygous in the solution
    /// * `skipped_variants` - the number of variants ignored in the solution
    /// # Panics
    /// * if `actual_cost < estimated_cost`, because that would imply something broken in our heuristics that makes A* no longer work
    pub fn astar_new(
        pruned_solutions: u64, estimated_cost: u64, actual_cost: u64,
        phased_variants: u64, phased_snvs: u64, homozygous_variants: u64, skipped_variants: u64
    ) -> PhaseStats {
        assert!(actual_cost >= estimated_cost);
        PhaseStats {
            pruned_solutions: Some(pruned_solutions),
            estimated_cost: Some(estimated_cost),
            actual_cost: Some(actual_cost),
            phased_variants: Some(phased_variants),
            phased_snvs: Some(phased_snvs),
            homozygous_variants: Some(homozygous_variants),
            skipped_variants: Some(skipped_variants)
        }
    }

    pub fn get_pruned_solutions(&self) -> Option<u64> {
        self.pruned_solutions
    }

    pub fn phased_snvs(&self) -> Option<u64> {
        self.phased_snvs
    }

    /// Returns the ratio of estimated_cost / actual_cost. In a perfect world, this is very near 1.0.
    pub fn get_cost_ratio(&self) -> Option<f64> {
        match self.estimated_cost {
            Some(ec) => {
                self.actual_cost.map(|ac| 
                    if ac == 0 {
                        assert_eq!(ec, 0);
                        1.0
                    } else {
                        ec as f64 / ac as f64
                    }
                )
            },
            None => None
        }
    }
}

/// This is a wrapper for writing out any stats to a file
pub struct StatsWriter {
    /// Handle for the CSV writer
    csv_writer: csv::Writer<File>
}

/// Contains all the data written to each row of our stats file
#[derive(Serialize)]
struct CsvRow {
    /// The index of the block
    block_index: usize,
    /// the sample for the block
    sample_name: String,
    /// the chromosome of the block
    chrom: String,
    /// the position of the first allele
    start: u64,
    /// the position of the last allele
    end: u64,
    /// the number of variants in the block
    num_variants: u64,
    /// the number of reads in the block
    num_reads: Option<u64>,
    /// the number of skipped reads in the block
    skipped_reads: Option<u64>,
    /// The number of variants loaded
    num_alleles: Option<u64>,
    /// Records the number of exact matches found for a type
    allele_matches: Option<String>,
    /// Records the number of inexact matches found for a type
    allele_partials: Option<String>,
    /// Records the number of no-matches found for a type
    allele_failures: Option<String>,
    /// Records the number of assignments to allele0 for a type
    allele0_assigned: Option<String>,
    /// Records the number of assignments to allele1 for a type
    allele1_assigned: Option<String>,
    /// The number of globally re-aligned reads
    global_aligned: Option<usize>,
    /// The number of locally re-aligned reads
    local_aligned: Option<usize>,
    /// The number of solutions that were pruned during calculation
    pruned_solutions: Option<u64>,
    /// For heuristic solvers, the estimate cost prior to computation
    estimated_cost: Option<u64>,
    /// The actual cost of the solution
    actual_cost: Option<u64>,
    /// The estimated / actual cost ratio
    cost_ratio: Option<f64>,
    /// The number of phased variants
    phased_variants: Option<u64>,
    /// The number of variants where the phasing solution turned them homozygous
    homozygous_variants: Option<u64>,
    /// The number of ignored variants
    skipped_variants: Option<u64>
}

impl StatsWriter {
    /// Creates a new writer for a given filename
    /// # Arguments
    /// * `filename` - the path to write all stats to
    pub fn new(filename: &Path) -> csv::Result<StatsWriter> {
        // modify the delimiter to "," if it ends with .csv
        let is_csv: bool = filename.extension().unwrap_or_default() == "csv";
        let delimiter: u8 = if is_csv { b',' } else { b'\t' };
        let csv_writer: csv::Writer<File> = csv::WriterBuilder::new()
            .delimiter(delimiter)
            .from_path(filename)?;
        Ok(StatsWriter {
            csv_writer
        })
    }

    /// Will write stats to a CSV file for us
    /// # Arguments
    /// * `phase_result` - the phasing results, which wraps block metadata, the read statistics, and the phasing statistics
    pub fn write_stats(&mut self, phase_result: &PhaseResult) -> csv::Result<()> {
        let num_reads;
        let skipped_reads;
        let num_alleles;
        let allele_matches;
        let allele_partials;
        let allele_failures;
        let allele0_assigned;
        let allele1_assigned;
        let global_aligned;
        let local_aligned;
        match &phase_result.read_statistics {
            Some(rs) => {
                num_reads = Some(rs.num_reads);
                skipped_reads = Some(rs.skipped_reads);
                num_alleles = Some(rs.num_alleles);
                allele_matches = Some(format!("{:?}", rs.exact_matches));
                allele_partials = Some(format!("{:?}", rs.inexact_matches));
                allele_failures = Some(format!("{:?}", rs.failed_matches));
                allele0_assigned = Some(format!("{:?}", rs.allele0_matches));
                allele1_assigned = Some(format!("{:?}", rs.allele1_matches));
                global_aligned = Some(rs.global_aligned);
                local_aligned = Some(rs.local_aligned);
            },
            None => {
                num_reads = None;
                skipped_reads = None;
                num_alleles = None;
                allele_matches = None;
                allele_partials = None;
                allele_failures = None;
                allele0_assigned = None;
                allele1_assigned = None;
                global_aligned = None;
                local_aligned = None;
            }
        };
        
        let pruned_solutions;
        let estimated_cost;
        let actual_cost;
        let cost_ratio;
        let phased_variants;
        let homozygous_variants;
        let skipped_variants;

        match &phase_result.statistics {
            Some(ps) => {
                pruned_solutions = ps.pruned_solutions;
                estimated_cost = ps.estimated_cost;
                actual_cost = ps.actual_cost;
                cost_ratio = ps.get_cost_ratio();
                phased_variants = ps.phased_variants;
                homozygous_variants = ps.homozygous_variants;
                skipped_variants = ps.skipped_variants;
            },
            None => {
                pruned_solutions = None;
                estimated_cost = None;
                actual_cost = None;
                cost_ratio = None;
                phased_variants = None;
                homozygous_variants = None;
                skipped_variants = None;
            }
        };

        let row: CsvRow = CsvRow {
            block_index: phase_result.phase_block.get_block_index(),
            sample_name: phase_result.phase_block.sample_name().to_string(),
            chrom: phase_result.phase_block.get_chrom().to_string(),
            start: phase_result.phase_block.get_start(),
            end: phase_result.phase_block.get_end(),
            num_variants: phase_result.phase_block.get_num_variants() as u64,
            num_reads,
            skipped_reads,
            num_alleles,
            allele_matches,
            allele_partials,
            allele_failures,
            allele0_assigned,
            allele1_assigned,
            global_aligned,
            local_aligned,
            pruned_solutions,
            estimated_cost,
            actual_cost,
            cost_ratio,
            phased_variants,
            homozygous_variants,
            skipped_variants
        };

        self.csv_writer.serialize(&row)?;
        self.csv_writer.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add_assign() {
        let mut read_stats: ReadStats = Default::default();
        let rhs = ReadStats::new(
            10, 12, 33, 
            [1; VariantType::Unknown as usize + 1],
            [2; VariantType::Unknown as usize + 1],
            [3; VariantType::Unknown as usize + 1],
            [2; VariantType::Unknown as usize + 1],
            [1; VariantType::Unknown as usize + 1],
            10, 2
        );

        println!("{read_stats:?}");
        read_stats += rhs.clone();
        println!("{rhs:?}");
        println!("{read_stats:?}");
        assert_eq!(read_stats, rhs);
    }
}
