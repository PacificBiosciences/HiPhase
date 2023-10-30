
use crate::phaser::HaplotagResult;

use log::{debug, trace, warn};
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use simple_error::bail;
use std::cell::RefCell;
use std::path::{PathBuf, Path};

/// Structure that maintains order of phase problems while writing solutions.
pub struct OrderedBamWriter {
    /// each bam writer is responsible for a specific sample
    sample_name: String,
    /// the template BAMs we are reading from and copying/modifying
    ref_bam_readers: Vec<RefCell<bam::IndexedReader>>,
    /// the BAMs we are writing to
    ref_bam_writers: Vec<RefCell<bam::Writer>>,
    /// the data that may be cached because we are waiting on earlier results
    map_store: HashMap<usize, HaplotagResult>,
    /// set of blocks that are indicated to skip
    skip_set: HashSet<usize>,
    /// the index of data we are waiting for
    current_index: usize,
    /// the most recently written chromosome
    current_chrom: String,
    /// the most recently written POS field
    current_pos: u64,
    /// tracks chromosomes that have been finalized
    finished_chroms: HashSet<String>
}

impl OrderedBamWriter {
    /// Creates a new `OrderedBamWriter` using template BAMs.
    /// # Arguments
    /// * `sample_name` - the sample this BAM writer corresponds to
    /// * `reference_filename` - the path to the reference genome
    /// * `input_bams` - the template VCF file containing unphased variants
    /// * `output_bams` - the VCF file that will get created containing our phase solutions
    /// * `thread_pool` - a shared thread pool for BAM I/O
    pub fn new(
        sample_name: String, reference_filename: &Path,
        input_bams: &[PathBuf], output_bams: &[PathBuf],
        thread_pool: &rust_htslib::tpool::ThreadPool
    ) -> Result<OrderedBamWriter, rust_htslib::errors::Error> {
        // log this
        debug!("Creating BAM writer for {}:", sample_name);
        debug!("\tInputs:  {:?}", input_bams);
        debug!("\tOutputs: {:?}", output_bams);

        // get all the stuff we need for reading the VCF setup
        let mut ref_bam_readers: Vec<RefCell<bam::IndexedReader>> = vec![];
        let mut ref_bam_writers: Vec<RefCell<bam::Writer>> = vec![];
        for (i, path) in input_bams.iter().enumerate() {
            let mut bam_reader: bam::IndexedReader = bam::IndexedReader::from_path(path)?;
            bam_reader.set_reference(reference_filename)?;
            bam_reader.set_thread_pool(thread_pool)?;
            let bam_header: bam::HeaderView = bam_reader.header().clone();
            let ref_bam_reader: RefCell<bam::IndexedReader> = RefCell::new(bam_reader);
            
            //now setup the outputs, we want to do the header stuff here also
            let mut output_header: bam::header::Header = bam::header::Header::from_template(&bam_header);
            let cli_string: String = std::env::args().collect::<Vec<String>>().join(" ");
            let cli_version: &str = &crate::cli::FULL_VERSION;

            let mut cli_record = bam::header::HeaderRecord::new("PG".as_bytes());
            cli_record.push_tag("PN".as_bytes(), &"hiphase");
            cli_record.push_tag("ID".as_bytes(), &format!("hiphase-v{cli_version}"));
            cli_record.push_tag("VN".as_bytes(), &cli_version);
            cli_record.push_tag("CL".as_bytes(), &cli_string);
            output_header.push_record(&cli_record);
            
            // TODO: do we need to add something to the BAM headers? I'm not sure yet, we'll come back to this
            // output_header.push_record(r#"##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">"#.as_bytes());
            let bam_format: bam::Format = if output_bams[i].extension().unwrap_or_default() == "cram" {
                bam::Format::Cram
            } else {
                bam::Format::Bam
            };
            let mut bam_writer: bam::Writer = bam::Writer::from_path(
                &output_bams[i],
                &output_header,
                bam_format
            )?;
            bam_writer.set_reference(reference_filename)?;
            bam_writer.set_thread_pool(thread_pool)?;
            let ref_bam_writer: RefCell<bam::Writer> = RefCell::new(bam_writer);

            // push everything into the vecs for storage
            ref_bam_readers.push(ref_bam_reader);
            ref_bam_writers.push(ref_bam_writer);
        }

        Ok(OrderedBamWriter {
            sample_name,
            ref_bam_readers,
            ref_bam_writers,
            map_store: Default::default(),
            skip_set: Default::default(),
            current_index: 0,
            current_chrom: "".to_string(),
            current_pos: 0,
            finished_chroms: Default::default()
        })
    }

    /// Returns the phase block result that the writer is currently waiting to receive.
    pub fn get_wait_block(&self) -> usize {
        self.current_index
    }

    /// Adds a phase block result to our queue for writing.
    /// # Arguments
    /// * `phase_result` - a phasing result that will be written in the correct order with other blocks
    pub fn write_phase_block(&mut self, haplotag_result: HaplotagResult) -> Result<(), Box<dyn std::error::Error>> {
        let block_index: usize = haplotag_result.phase_block.get_block_index();
        if block_index < self.current_index {
            bail!("Block index is smaller than next expected index");
        }
        if haplotag_result.phase_block.sample_name() != self.sample_name {
            bail!("Received haplotag result for sample other than the one specified");
        }
        match self.map_store.insert(block_index, haplotag_result) {
            None => {},
            Some(_) => {
                bail!("Block index was already present in the map_store");
            }
        };
        self.drain_map_store()?;
        Ok(())
    }

    /// Adds a dummy block result to our queue for writing.
    /// This is only necessary because we have block indices that are shared across samples, but this only cares about one sample.
    /// # Arguments
    /// * `phase_result` - a phasing result that will be written in the correct order with other blocks
    pub fn write_dummy_block(&mut self, block_index: usize) -> Result<(), Box<dyn std::error::Error>> {
        if block_index < self.current_index {
            bail!("Block index is smaller than next expected index");
        }

        // add to the skip list and drain away
        self.skip_set.insert(block_index);
        self.drain_map_store()?;
        Ok(())
    }

    /// This will drain phase block solutions in the correct order if they have been received.
    /// It drains as far as it can given the current results and then stops to wait for more data.
    fn drain_map_store(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        while !self.map_store.is_empty() {
            match self.map_store.remove(&self.current_index) {
                Some(haplotag_result) => {
                    trace!("Draining {}", self.current_index);

                    let chrom_result: &str = haplotag_result.phase_block.get_chrom();
                    if chrom_result != self.current_chrom {
                        if self.current_chrom.is_empty() {
                            // this is the first block, so lets set chrom and move on
                            self.current_chrom = chrom_result.to_string();
                        } else {
                            // the next block is on a different chromosome, so we need to finalize this chromosome
                            self.finalize_chromosome()?;

                            // now setup for the next chromosome
                            self.current_chrom = chrom_result.to_string();
                            self.current_pos = 0;
                        }
                    } else {
                        // the chromosome matches previous, I don't think we actually have to do anything
                        // as long as our assumptions are correct
                    }

                    /*
                    General plan: loop over the writers and write each BAM separately
                    - get all of the index-specific readers/writers ready
                    - adjust haplotype_index to only correspond to the current vcf_index
                    - we may need to adjust our errors to allow for empty blocks; we can make a dummy VCF to test with
                    - ^ we should make a dummy VCF anyways just to verify that nothing in our output changes by adding an empty VCF
                    */
                    let start_pos = self.current_pos;
                    let end_pos = haplotag_result.phase_block.get_end();
                        
                    for (bam_index, ref_bam_writer) in self.ref_bam_writers.iter().enumerate() {
                        // prep the writer
                        let mut bam_writer = ref_bam_writer.borrow_mut();
                        
                        // we _should_ just have to iterate from the last written through the end of the block
                        let mut bam_reader = self.ref_bam_readers[bam_index].borrow_mut();
                        
                        // key is a read name, value is (block_id, haplotype)
                        let read_block_lookup = &haplotag_result.reads;
                        
                        // annoyingly, the bam reader (here) uses normal rust conventions for 0..len (end is exclusive), whereas the VCF reader is inclusive
                        // so to fix, we need to add the +1 as shown below
                        match bam_reader.fetch(bam::FetchDefinition::RegionString(chrom_result.as_bytes(), start_pos as i64, end_pos as i64 + 1)) {
                            Ok(()) => {
                                for record_result in bam_reader.records() {
                                    // set up the record for the new BAM file
                                    let mut record = record_result?;
                                    let record_pos = record.pos();
                                    if record_pos < start_pos as i64 {
                                        // this can happen when you have reads that overlap the location but don't start *after* the position
                                        // we have already written though, so don't write it again
                                        continue;
                                    }

                                    // this may need to be <=, hard to tell yet
                                    assert!(record_pos <= end_pos as i64);

                                    // now check if the read name has a lookup
                                    let read_name = std::str::from_utf8(record.qname()).unwrap();
                                    match read_block_lookup.get(read_name) {
                                        Some((phase_block_id, haplotag)) => {
                                            // we have a match, modify phase info
                                            // phase_block_id is 0-based, so add 1 to it
                                            record.push_aux("PS".as_bytes(), bam::record::Aux::I32((phase_block_id + 1).try_into()?))?;
                                            // haplotag is 0/1 and we want 1/2 in the BAM, so add 1 to it
                                            record.push_aux("HP".as_bytes(), bam::record::Aux::U8((haplotag + 1).try_into()?))?;
                                            bam_writer.write(&record)?;
                                        },
                                        None => {
                                            // no match, so just copy the read over
                                            bam_writer.write(&record)?;
                                        }
                                    };
                                }
                            },
                            Err(e) => {
                                if end_pos == 0 {
                                    warn!("Empty problem block received, no read mappings on chromosome {}", chrom_result);
                                } else {
                                    warn!("Received \'{}\', while seeking to {}:{}-{} in bam #{}, likely no reads in region", e, chrom_result, start_pos, end_pos, bam_index);
                                    //return Err(e);
                                }
                            }
                        };
                    }

                    // set up for the next block we get
                    self.current_pos = end_pos+1;
                    self.current_index += 1;
                },
                None => {
                    // it's not in the received blocks, check if it is in the skip set
                    if self.skip_set.remove(&self.current_index) {
                        // we did find it here
                        self.current_index += 1;
                    } else {
                        // it's not there either, time to end the looping for now
                        break;
                    }
                }
            };
        }
        Ok(())
    }

    /// Writes all remaining variants for a chromosome.
    /// It should be called after all phase blocks for a chromosome are received.
    /// # Panics
    /// * if the current chromosome has been previously finalized
    pub fn finalize_chromosome(&mut self) -> Result<(), rust_htslib::errors::Error> {
        // make sure we haven't done this chromosome before
        assert!(!self.finished_chroms.contains(&self.current_chrom));

        // finalize the area
        let start_pos = self.current_pos;
        for (bam_index, ref_bam_writer) in self.ref_bam_writers.iter().enumerate() {
            // prep the writer
            let mut bam_writer = ref_bam_writer.borrow_mut();
            let mut bam_reader = self.ref_bam_readers[bam_index].borrow_mut();
            
            //match vcf_reader.fetch(chrom_index, start_pos, None) {
            match bam_reader.fetch(bam::FetchDefinition::RegionString(self.current_chrom.as_bytes(), start_pos as i64, i64::MAX)) {
                Ok(()) => {
                    for record_result in bam_reader.records() {
                        // set up the record for the new BAM file
                        let record = record_result?;
                        let record_pos = record.pos();
                        if record_pos < start_pos as i64 {
                            // this can happen when you have reads that overlap the location but don't start *after* the position
                            // we have already written though, so don't write it again
                            continue;
                        }

                        // nothing left should be tagged
                        bam_writer.write(&record)?;
    
                        // adding this last bit should prevent double writes by accident from a user
                        self.current_pos = self.current_pos.max(record_pos as u64 + 1);
                    }
                },
                Err(e) => {
                    warn!("Received \'{}\', likely caused by no trailing variants detected for {}:{}-END", e, self.current_chrom, start_pos);
                }
            };
        }

        self.finished_chroms.insert(self.current_chrom.clone());
        Ok(())
    }

    pub fn copy_remaining_chromosomes(&mut self) -> Result<(), rust_htslib::errors::Error> {
        // go through each BAM and just copy any chromosomes we didn't do
        for (bam_index, ref_bam_writer) in self.ref_bam_writers.iter().enumerate() {
            // prep the writer
            let mut bam_writer = ref_bam_writer.borrow_mut();
            let mut bam_reader = self.ref_bam_readers[bam_index].borrow_mut();
            debug!("Finalizing bam #{}...", bam_index);

            // first, find anything that didn't copy
            let mut remaining_contigs: Vec<bam::FetchDefinition> = vec![];
            let mut remaining_contigs_copy: Vec<bam::FetchDefinition> = vec![]; // we have to do this because there is no clone for FetchDefs
            let target_names: Vec<String> = bam_reader.header().target_names().iter()
                .map(|s| String::from_utf8(s.to_vec()).unwrap())
                .collect();
            for tn in target_names.iter() {
                if !self.finished_chroms.contains(tn) {
                    // this target has not been completed yet
                    remaining_contigs.push(bam::FetchDefinition::String(tn.as_bytes()));
                    remaining_contigs_copy.push(bam::FetchDefinition::String(tn.as_bytes()));
                }
            }

            // add unmapped at the end
            remaining_contigs.push(bam::FetchDefinition::Unmapped);
            remaining_contigs_copy.push(bam::FetchDefinition::Unmapped);

            debug!("Remaining contigs: {:?}", remaining_contigs);

            // now go through each, and fully copy everything
            for (fetch_index, fetch_definition) in remaining_contigs.into_iter().enumerate() {
                match bam_reader.fetch(fetch_definition) {
                    Ok(()) => {
                        for record_result in bam_reader.records() {
                            // set up the record for the new BAM file
                            let record = record_result?;
                            
                            // nothing left should be tagged
                            bam_writer.write(&record)?;
                        }
                    },
                    Err(e) => {
                        warn!("Received \'{}\', likely caused by no reads detected for {:?}", e, remaining_contigs_copy[fetch_index]);
                    }
                };
            }
        }

        // if we make it here, everything should be good
        Ok(())
    }
}
