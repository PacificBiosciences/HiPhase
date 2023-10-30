
use serde::Serialize;
use std::fs::File;
use std::path::Path;

use crate::phaser::HaplotagResult;

/// This is a wrapper for writing out any stats to a file
pub struct HaplotagWriter {
    /// Handle for the CSV writer
    csv_writer: csv::Writer<File>
}

/// Contains all the data written to each row of our stats file
#[derive(Serialize)]
struct HaplotagRow {
    /// The index of the block
    source_block_index: usize,
    /// the sample name tied to the block
    sample_name: String,
    /// the chromosome of the block
    chrom: String,
    /// Phase set ID - usually the position first variant in the block
    phase_block_id: u64,
    /// the read name that is assigned
    read_name: String,
    /// the haplotype the read is assigned to
    haplotag: u8
}

impl HaplotagWriter {
    /// Creates a new writer for a given filename
    /// # Arguments
    /// * `filename` - the path to write all stats to
    pub fn new(filename: &Path) -> csv::Result<HaplotagWriter> {
        // modify the delimiter to "," if it ends with .csv
        let is_csv: bool = filename.extension().unwrap_or_default() == "csv";
        let delimiter: u8 = if is_csv { b',' } else { b'\t' };
        let csv_writer: csv::Writer<File> = csv::WriterBuilder::new()
            .delimiter(delimiter)
            .from_path(filename)?;
        Ok(HaplotagWriter {
            csv_writer
        })
    }

    /// Writes the haplotag results from a given block.
    /// # Arguments
    /// * `haplotag_result` - the HaplotagResult from phasing, in-order is not required
    /// # Errors
    /// * if the csv_writer has any errors
    pub fn write_block(&mut self, haplotag_result: &HaplotagResult) -> Result<(), Box<dyn std::error::Error>> {
        let source_block_index: usize = haplotag_result.phase_block.get_block_index();
        let sample_name: String = haplotag_result.phase_block.sample_name().to_string();
        let chrom: String = haplotag_result.phase_block.get_chrom().to_string();
        for (read_name, &(phase_block_id_0, haplotag_0)) in haplotag_result.reads.iter() {
            let phase_block_id: u64 = (phase_block_id_0 + 1).try_into()?;
            let haplotag: u8 = (haplotag_0 + 1).try_into()?;
            let row: HaplotagRow = HaplotagRow {
                source_block_index,
                sample_name: sample_name.clone(),
                chrom: chrom.clone(),
                phase_block_id,
                read_name: read_name.clone(),
                haplotag
            };
            self.csv_writer.serialize(&row)?;
            self.csv_writer.flush()?;
        }
        Ok(())
    }
}