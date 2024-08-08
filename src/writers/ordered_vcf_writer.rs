
use crate::block_gen::is_phasable_variant;
use crate::phaser::PhaseResult;

use log::{debug, trace};
use rust_htslib::bcf;
use rust_htslib::bcf::Read;
use rust_htslib::bcf::record::GenotypeAllele;
use rustc_hash::FxHashMap as HashMap;
use simple_error::bail;
use std::cell::RefCell;
use std::collections::VecDeque;
use std::io;
use std::path::PathBuf;

/// Mostly for keeping the necessary information for a phase together when we store the results in memory prior to writing
#[derive(Debug)]
struct SingleVariantPhase {
    /// haplotype 1
    h1: u8,
    /// haplotype 2
    h2: u8,
    /// phase block ID
    block_id: usize
}

/// Structure that maintains order of phase problems while writing solutions.
pub struct OrderedVcfWriter {
    /// the template VCF we are reading from and copying if no phases are given
    ref_vcf_readers: Vec<RefCell<bcf::IndexedReader>>,
    /// A copy of the input VCF header, cached here for performance
    vcf_headers: Vec<bcf::header::HeaderView>,
    /// the VCF we are writing to
    ref_vcf_writers: Vec<RefCell<bcf::Writer>>,
    /// the data that may be cached because we are waiting on earlier results
    map_store: HashMap<usize, PhaseResult>,
    /// the index of data we are waiting for
    current_index: usize,
    /// the most recently written chromosome
    current_chrom: String,
    /// the most recently written POS field
    current_pos: u64,
    /// the distance you can write for each sample
    current_positions: HashMap<String, u64>,
    /// The minimum quality for variants that were phased. If this is different from what is generating blocks, there will likely be a panic.
    min_quality: i32,
    /// Per VCF file, there is a hashmap of sample_name -> index in VCF file
    sample_indices: Vec<HashMap<String, usize>>,
    /// Per VCF file, there is a hashmap of sample_name -> queue of variants
    phase_queues: Vec<HashMap<String, VecDeque<SingleVariantPhase>>>
}

impl OrderedVcfWriter {
    /// Creates a new `OrderedVcfWriter` using template VCFs.
    /// # Arguments
    /// * `input_vcfs` - the template VCF files containing unphased variants
    /// * `output_vcfs` - the VCF files that will get created containing our phase solutions
    /// * `min_quality` - the minimum quality that indicates a variant to phase
    /// * `sample_name` - the sample getting phased
    pub fn new(input_vcfs: &[PathBuf], output_vcfs: &[PathBuf], min_quality: i32, sample_names: &[String]) -> Result<OrderedVcfWriter, Box<dyn std::error::Error>> {
        //get all the stuff we need for reading the VCF setup
        let mut ref_vcf_readers: Vec<RefCell<bcf::IndexedReader>> = vec![];
        let mut vcf_headers: Vec<bcf::header::HeaderView> = vec![];
        let mut ref_vcf_writers: Vec<RefCell<bcf::Writer>> = vec![];
        let mut sample_indices: Vec<HashMap<String, usize>> = vec![];
        let mut phase_queues: Vec<HashMap<String, VecDeque<SingleVariantPhase>>> = vec![];
        for (i, path) in input_vcfs.iter().enumerate() {
            let vcf_reader: bcf::IndexedReader = bcf::IndexedReader::from_path(path)?;
            let vcf_header: bcf::header::HeaderView = vcf_reader.header().clone();
            let ref_vcf_reader: RefCell<bcf::IndexedReader> = RefCell::new(vcf_reader);

            // first make sure we find the sample in this file
            let mut vcf_sample_hash: HashMap<String, usize> = Default::default();
            let mut vcf_phase_queue: HashMap<String, VecDeque<SingleVariantPhase>> = Default::default();
            for sample_name in sample_names.iter() {
                let mut lookup_index: Option<usize> = None;
                for (sample_index, &vcf_sample) in vcf_header.samples().iter().enumerate() {
                    let vcf_sample_string: String = std::str::from_utf8(vcf_sample).unwrap().to_string();
                    if &vcf_sample_string == sample_name {
                        lookup_index = Some(sample_index);
                        break;
                    }
                }
                let lookup_index: usize = match lookup_index {
                    Some(index) => {
                        index
                    },
                    None => {
                        bail!("Sample name {:?} was not found in VCF: {:?}", sample_name, path);
                    }
                };

                // add the sample lookup and initialize to an empty queue
                vcf_sample_hash.insert(sample_name.to_string(), lookup_index);
                vcf_phase_queue.insert(sample_name.to_string(), Default::default());
            }
            
            //now setup the outputs, we want to do the header stuff here also
            let mut output_header: bcf::header::Header = bcf::header::Header::from_template(&vcf_header);
            let cli_string: String = std::env::args().collect::<Vec<String>>().join(" ");
            let cli_version: &str = &crate::cli::FULL_VERSION;
            output_header.push_record(format!(r#"##HiPhase_version="{cli_version}""#).as_bytes());
            output_header.push_record(format!(r#"##HiPhase_command="{cli_string}""#).as_bytes());
            output_header.push_record(r#"##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">"#.as_bytes());
            output_header.push_record(r#"##FORMAT=<ID=PF,Number=1,Type=String,Description="Phasing flag">"#.as_bytes());
            let vcf_writer: bcf::Writer = bcf::Writer::from_path(
                &output_vcfs[i],
                &output_header,
                false,
                bcf::Format::Vcf
            )?;
            let ref_vcf_writer: RefCell<bcf::Writer> = RefCell::new(vcf_writer);

            // push everything into the vecs for storage
            ref_vcf_readers.push(ref_vcf_reader);
            vcf_headers.push(vcf_header);
            ref_vcf_writers.push(ref_vcf_writer);
            sample_indices.push(vcf_sample_hash);
            phase_queues.push(vcf_phase_queue);
        }

        // default all of them to start at 0
        let mut current_positions: HashMap<String, u64> = Default::default();
        for sample in sample_names.iter() {
            current_positions.insert(sample.clone(), 0);
        }

        Ok(OrderedVcfWriter {
            ref_vcf_readers,
            vcf_headers,
            ref_vcf_writers,
            map_store: Default::default(),
            current_index: 0,
            current_chrom: "".to_string(),
            current_pos: 0,
            current_positions,
            min_quality,
            sample_indices,
            phase_queues
        })
    }

    /// Returns the phase block result that the writer is currently waiting to receive.
    pub fn get_wait_block(&self) -> usize {
        self.current_index
    }

    /// Adds a phase block result to our queue for writing.
    /// # Arguments
    /// * `phase_result` - a phasing result that will be written in the correct order with other blocks
    pub fn write_phase_block(&mut self, phase_result: PhaseResult) -> Result<(), Box<dyn std::error::Error>> {
        let block_index: usize = phase_result.phase_block.get_block_index();
        if block_index < self.current_index {
            return Err(Box::new(io::Error::new(io::ErrorKind::Other, "Block index is smaller than next expected index")));
        }
        match self.map_store.insert(block_index, phase_result) {
            None => {},
            Some(_) => {
                return Err(Box::new(io::Error::new(io::ErrorKind::Other, "Block index was already present in the map_store")));
            }
        };
        self.drain_map_store()
    }

    /// This will drain phase block solutions in the correct order if they have been received.
    /// It drains as far as it can given the current results and then stops to wait for more data.
    fn drain_map_store(&mut self) -> Result<(), Box<dyn std::error::Error>> {

        /*
         * New plan for multi-sample VCFs:
         * 1. we check if we have the next phase result like we already do
         * 2. we append all of the relevant info from that phase block into the variant queue for a given sample
         *   - do we have one queue per sample? or one per sample-VCF combination?
         *   - probably sample-VCF combination
         * 3. we iterate over each VCF from current position to the minimum end position of all sample blocks so far
         *   - this _should_ guarantee at least one sample queue is empty
         *   - for each variant, we check if it's getting phase for the sample and modify if so
         *   - this will pop an entry from the variant queue for the sample
         * 4. at the end of a chromosome:
         *   - iterate as normal
         *   - check that all queues are empty
         *   - put the new info into the queue, rinse and repeat
         */

        while !self.map_store.is_empty() {
            match self.map_store.remove(&self.current_index) {
                Some(phase_result) => {
                    trace!("Draining {}", self.current_index);

                    let chrom_result: &str = phase_result.phase_block.get_chrom();
                    if chrom_result != self.current_chrom {
                        if self.current_index == 0 {
                            // this is the first block, so lets set chrom and move on
                            self.current_chrom = chrom_result.to_string();
                        } else {
                            // the next block is on a different chromosome, so we need to finalize this chromosome
                            self.write_to_end_position()?;

                            // now setup for the next chromosome
                            self.current_chrom = chrom_result.to_string();
                            self.current_pos = 0;
                            for (_k, v) in self.current_positions.iter_mut() {
                                // set all current positions back to 0
                                *v = 0;
                            }
                        }
                    } else {
                        // the chromosome matches previous, I don't think we actually have to do anything
                        // as long as our assumptions are correct
                    }

                    let sample_name = phase_result.phase_block.sample_name();
                    for (vcf_index, phase_queue) in self.phase_queues.iter_mut().enumerate() {
                        let sample_queue: &mut VecDeque<SingleVariantPhase> = phase_queue.get_mut(sample_name).unwrap();
                        let mut previous_block_id: usize = 0;
                        for (haplotype_index, &h1_index) in phase_result.haplotype_1.iter().enumerate() {
                            if vcf_index == phase_result.variants[haplotype_index].get_vcf_index() {
                                // h1 and h2 are just internal representations
                                let h2_index: u8 = phase_result.haplotype_2[haplotype_index];

                                // convert them to the file representations where possible
                                let h1 = phase_result.variants[haplotype_index].convert_index(h1_index);
                                let h2 = phase_result.variants[haplotype_index].convert_index(h2_index);

                                // add one here because we need it to be 1-based
                                let block_id: usize = phase_result.block_ids[haplotype_index]+1;
                                if haplotype_index == 0 || block_id != previous_block_id {
                                    debug!("New block ID found for {}: {}", self.current_index, block_id);
                                }
                                previous_block_id = block_id;

                                sample_queue.push_back(SingleVariantPhase { h1, h2, block_id });
                            } else {
                                // this variant is not a part of this VCF file
                            }
                        }
                    }

                    // update the minimum position for this sample and then tell the writer to do what it can
                    *self.current_positions.get_mut(sample_name).unwrap() = phase_result.phase_block.get_end();
                    self.write_to_min_position()?;
                    self.current_index += 1;
                },
                None => {
                    break;
                }
            };
        }
        Ok(())
    }

    /// This will trigger the writer to try to write everything that remains on this chromosome
    /// # Errors
    /// * if we have any issues writing the VCF file
    /// * if when we finish writing, there are still blocks in the queue
    pub fn write_to_end_position(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        self.write_to_position(u64::MAX)?;
        
        // sanity check to make sure all queues are empty once we get to this point
        for hm in self.phase_queues.iter() {
            for (_sample_name, queue) in hm.iter() {
                if !queue.is_empty() {
                    bail!("Finished writing chromosome, but variant queues are not empty");
                }
            }
        }

        Ok(())
    }

    /// This will trigger to writer to write everything up to the minimum position on this chromosome
    /// # Errors
    /// * if we have any issues writing the VCF file
    fn write_to_min_position(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        let min_position: u64 = *self.current_positions.values().min().unwrap();
        self.write_to_position(min_position)
    }

    /// This should write out all VCF lines up to a certain position and then break.
    /// # Arguments
    /// * `final_position` - the last position that is included in this write
    /// # Errors
    /// * if we have any issues writing the VCF file
    fn write_to_position(&mut self, final_position: u64) -> Result<(), Box<dyn std::error::Error>> {
        if self.current_pos == final_position {
            // we can't write anything new
            return Ok(());
        }

        // 
        let start_pos = self.current_pos;
        for (vcf_index, ref_vcf_writer) in self.ref_vcf_writers.iter().enumerate() {
            // prep the writer
            let mut vcf_writer = ref_vcf_writer.borrow_mut();
            let mut vcf_reader = self.ref_vcf_readers[vcf_index].borrow_mut();
            let chrom_index: u32 = self.vcf_headers[vcf_index].name2rid(self.current_chrom.as_bytes())?;
            // let sample_index: usize = self.sample_indices[vcf_index];
            
            match vcf_reader.fetch(chrom_index, start_pos, Some(final_position)) {
                Ok(()) => {
                    for record_result in vcf_reader.records() {
                        // set up the record for the new VCF file
                        let mut record = record_result?;
                        vcf_writer.translate(&mut record);
                        let record_pos = record.pos();
                        if record_pos < start_pos as i64 {
                            // this can happen when you have very very long indels that span one of our breaks
                            // we have already written though, so don't write it again
                            continue;
                        }
                        
                        // we now have to iterate over each sample and modify the entries as necessary
                        let vcf_sample_indices = &self.sample_indices[vcf_index];
                        let mut changes_made: bool = false;

                        // initialize the alleles array to be the same, these may change inside the loop
                        let mut alleles = vec![];
                        let record_gt = record.genotypes().unwrap();
                        for si in 0..record.sample_count() {
                            let genotype = record_gt.get(si as usize);
                            match genotype.len() {
                                0 => bail!("Encountered empty genotype record at position {}", record.pos()),
                                1 => {
                                    // TRGT can make single-allele GT calls, just copy it over as normal
                                    // it will not be modified below because it is a homozygous allele
                                    alleles.push(i32::from(genotype[0]));
                                    alleles.push(i32::MIN+1); // sentinel value for end
                                },
                                2 => {
                                    // this is 99.99999999% path
                                    alleles.push(i32::from(genotype[0]));
                                    alleles.push(i32::from(genotype[1]));
                                },
                                gt_len => {
                                    // we do not have 3+ GT entries implemented
                                    bail!("Encountered GT of length {} at record {}", gt_len, record.desc())
                                }
                            }
                        }

                        // initially empty PS blocks also
                        let mut ps_blocks: Vec<Vec<u8>> = vec![".".as_bytes().to_vec(); record.sample_count() as usize];
                        let mut phase_flags: Vec<Vec<u8>> = vec![".".as_bytes().to_vec(); record.sample_count() as usize];
                        let mut flagged_variants: bool = false;

                        for (sample_name, &sample_index) in vcf_sample_indices.iter() {
                            // now see how we handle this variant
                            let include_variant = is_phasable_variant(&record, sample_index, self.min_quality, false)?;
                            if include_variant {
                                // sanity checks
                                assert!(u64::try_from(record_pos).unwrap() >= start_pos);
                                assert!(u64::try_from(record_pos).unwrap() <= final_position);

                                let variant_to_write = match self.phase_queues[vcf_index].get_mut(sample_name).unwrap().pop_front() {
                                    Some(v) => v,
                                    None => {
                                        bail!("Variant requested from empty queue during VCF writing");
                                    }
                                };

                                // these are already converted to the VCF entries, so just compare to see if we need to overwrite
                                let h1 = variant_to_write.h1;
                                let h2 = variant_to_write.h2;
                                if h1 == h2 {
                                    // algorithm decided it was better if these were homozygous allele
                                    // for now, we will just write out the original record

                                    if h1 == u8::MAX {
                                        // these were intentionally ignored by HiPhase, mark it as such
                                        phase_flags[sample_index] = "TR_OVERLAP".as_bytes().to_vec();
                                        flagged_variants = true;
                                    }

                                } else {
                                    // we need to alter the genotypes for this sample to phased
                                    let sample_gt_offset: usize = 2 * sample_index;
                                    alleles[sample_gt_offset] = i32::from(GenotypeAllele::Unphased(h1 as i32));
                                    alleles[sample_gt_offset + 1] = i32::from(GenotypeAllele::Phased(h2 as i32));

                                    // the push_format_string expects &[u8] bytes so we have to:
                                    //   1. convert the output to a String
                                    //   2. interpret that to bytes
                                    //   3. convert to a Vec for ownership
                                    ps_blocks[sample_index] = variant_to_write.block_id
                                        .to_string().as_bytes().to_vec();
                                    changes_made = true;
                                }
                            } else {
                                // this variant is not included in phasing, so we can just leave it as is
                            }
                        }

                        if changes_made {
                            // if we altered something, then alter the record and add PS
                            record.push_format_integer(b"GT", &alleles)?;
                            record.push_format_string("PS".as_bytes(), &ps_blocks).unwrap();
                        }
                        if flagged_variants {
                            // we have at least one variant that was ignored, use PF = PhaseFlag 
                            record.push_format_string("PF".as_bytes(), &phase_flags).unwrap();
                        }
                        // all modifications have been made, write it out
                        vcf_writer.write(&record)?;
                    }
                },
                Err(e) => {
                    if final_position == 0 {
                        debug!("Empty problem block received, no heterozygous variants on chromosome {}", self.current_chrom);
                    } else {
                        debug!("Received \'{}\', while seeking to {}:{}-{} in vcf #{}, likely no variants present", e, self.current_chrom, start_pos, final_position, vcf_index);
                    }
                }
            }
        }
        
        // we wrote out things at final_position, so go one past it
        if final_position == u64::MAX {
            self.current_pos = final_position;
        } else {
            self.current_pos = final_position+1;
        }
        Ok(())
    }
}
