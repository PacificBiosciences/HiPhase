
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
            append_hiphase_program_record(&mut output_header, &bam_header, cli_version, &cli_string);
            
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

                                    // quick sanity check
                                    assert!(record_pos <= end_pos as i64);

                                    // no matter what, we need to strip any existing phase information from the record
                                    strip_record_phasing(&mut record)?;

                                    // now check if the read name has a lookup
                                    let read_name = std::str::from_utf8(record.qname()).unwrap();

                                    if let Some((phase_block_id, haplotag)) = read_block_lookup.get(read_name) {
                                        // we have a match, modify phase info
                                        // phase_block_id is 0-based, so add 1 to it
                                        record.push_aux("PS".as_bytes(), bam::record::Aux::I32((phase_block_id + 1).try_into()?))?;
                                        // haplotag is 0/1 and we want 1/2 in the BAM, so add 1 to it
                                        record.push_aux("HP".as_bytes(), bam::record::Aux::U8((haplotag + 1).try_into()?))?;
                                    } else {
                                        // no haplotag information for this read, the record was already stripped; no-op
                                    };
                                    bam_writer.write(&record)?;
                                }
                            },
                            Err(e) => {
                                if end_pos == 0 {
                                    warn!("Empty problem block received, no read mappings on chromosome {}", chrom_result);
                                } else {
                                    warn!("Received \'{}\', while seeking to {}:{}-{} in bam #{}, likely no reads in region", e, chrom_result, start_pos, end_pos, bam_index);
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
                        let mut record = record_result?;
                        let record_pos = record.pos();
                        if record_pos < start_pos as i64 {
                            // this can happen when you have reads that overlap the location but don't start *after* the position
                            // we have already written though, so don't write it again
                            continue;
                        }

                        // nothing left should be tagged
                        strip_record_phasing(&mut record)?;
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
                            let mut record = record_result?;
                            
                            // nothing left should be tagged
                            strip_record_phasing(&mut record)?;
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

/// Collects existing `@PG` ID values from a BAM header.
/// Parses `@PG` lines from the header text so comment (`@CO`) lines are ignored.
/// # Arguments
/// * `header` - the BAM header whose text is scanned for `@PG` ID fields
fn collect_pg_ids(header: &bam::HeaderView) -> HashSet<String> {
    // TODO: this method converts the header to bytes and then parses like a string
    //       it would be better to parse the header directly, but machinery is missing from htslib
    let mut ids: HashSet<String> = Default::default();
    for line in header.as_bytes().split(|&b| b == b'\n') {
        if !line.starts_with(b"@PG") {
            continue;
        }
        for field in line.split(|&b| b == b'\t') {
            if let Some(id_bytes) = field.strip_prefix(b"ID:") {
                if let Ok(id) = std::str::from_utf8(id_bytes) {
                    if !id.is_empty() {
                        ids.insert(id.to_string());
                    }
                }
            }
        }
    }
    ids
}

/// Returns `base_id`, or `base_id.1`, `base_id.2`, ... if that ID is already present.
/// # Arguments
/// * `existing_ids` - `@PG` IDs already present in the header
/// * `base_id` - the preferred program ID for this HiPhase run
fn unique_pg_id(existing_ids: &HashSet<String>, base_id: &str) -> String {
    if !existing_ids.contains(base_id) {
        return base_id.to_string();
    }
    let mut suffix: u32 = 1;
    loop {
        let candidate = format!("{base_id}.{suffix}");
        if !existing_ids.contains(&candidate) {
            return candidate;
        }
        suffix += 1;
    }
}

/// Appends a HiPhase `@PG` record with an ID that does not collide with records already in the header.
/// # Arguments
/// * `output_header` - the BAM header that will receive the new `@PG` record
/// * `template_header` - the input BAM header whose existing `@PG` records are scanned for ID collisions
/// * `cli_version` - the HiPhase version string written to the `VN` tag and used in the program ID
/// * `cli_string` - the command line written to the `CL` tag
fn append_hiphase_program_record(
    output_header: &mut bam::header::Header,
    template_header: &bam::HeaderView,
    cli_version: &str,
    cli_string: &str,
) {
    // figure out a unique program ID that doesn't collide with existing ones
    let base_id = format!("hiphase-v{cli_version}");
    let existing_ids = collect_pg_ids(template_header);
    let pg_id = unique_pg_id(&existing_ids, &base_id);

    // build the new record
    let mut cli_record = bam::header::HeaderRecord::new(b"PG");
    cli_record.push_tag(b"PN", &"hiphase");
    cli_record.push_tag(b"ID", &pg_id);
    cli_record.push_tag(b"VN", &cli_version);
    cli_record.push_tag(b"CL", &cli_string);
    output_header.push_record(&cli_record);
}

/// Removes the phase tags that HiPhase will normally write from a record.
/// This is so we erase anything that did not come from the tool.
fn strip_record_phasing(record: &mut bam::Record) -> Result<(), rust_htslib::errors::Error> {
    // much simpler than the VCF version since there is no multi-sample to worry about, just drop these tags
    match record.remove_aux(b"HP") {
        // either it was remove or was not there to begin with
        Ok(()) |
        Err(rust_htslib::errors::Error::BamAuxTagNotFound) => {},
        // some other crazy error
        Err(e) => {
            return Err(e);
        }
    };
    match record.remove_aux(b"PS") {
        // either it was remove or was not there to begin with
        Ok(()) |
        Err(rust_htslib::errors::Error::BamAuxTagNotFound) => Ok(()),
        // some other crazy error
        Err(e) => Err(e)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::HeaderView;
    
    /// Simple empty bam record creator; if we just use a "new()" record, htslib is unhappy when we access Aux
    fn dummy_bam_record() -> bam::Record {
        let mut record = bam::Record::new();
        record.set(
            b"test_qname",
            None,
            b"A",
            &[255]
        );
        record
    }

    /// Builds a minimal SAM header with optional `@PG` lines for uniqueness tests.
    fn header_with_pg_ids(pg_ids: &[&str]) -> bam::header::Header {
        let mut header = bam::header::Header::new();
        let mut hd = bam::header::HeaderRecord::new(b"HD");
        hd.push_tag(b"VN", &"1.6");
        hd.push_tag(b"SO", &"coordinate");
        header.push_record(&hd);

        let mut sq = bam::header::HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", &"chr1");
        sq.push_tag(b"LN", &"1000");
        header.push_record(&sq);

        for pg_id in pg_ids {
            let mut pg = bam::header::HeaderRecord::new(b"PG");
            pg.push_tag(b"ID", pg_id);
            pg.push_tag(b"PN", &"hiphase");
            header.push_record(&pg);
        }
        header
    }

    /// Appends a HiPhase `@PG` record and returns the header after an htslib round-trip.
    fn append_and_roundtrip(prior_pg_ids: &[&str], version: &str) -> HeaderView {
        let template = header_with_pg_ids(prior_pg_ids);
        let template_view = HeaderView::from_header(&template);
        let mut output_header = bam::header::Header::from_template(&template_view);
        append_hiphase_program_record(&mut output_header, &template_view, version, "hiphase --test");
        HeaderView::from_header(&output_header)
    }

    fn pg_ids_in_header(header: &HeaderView) -> Vec<String> {
        collect_pg_ids(header).into_iter().collect()
    }

    #[test]
    fn test_clear_record() {
        // make sure running on a record without the tags is fine (this will be *most* cases)
        let mut record = dummy_bam_record();
        strip_record_phasing(&mut record).unwrap();

        // create a record with the tags
        let mut record = dummy_bam_record();
        record.push_aux(b"HP", bam::record::Aux::U8(1)).unwrap();
        record.push_aux(b"PS", bam::record::Aux::I32(42)).unwrap();
        assert!(record.aux(b"HP").is_ok());
        assert!(record.aux(b"PS").is_ok());

        // strip them and verify they are gone
        strip_record_phasing(&mut record).unwrap();
        assert_eq!(record.aux(b"HP").err().unwrap(), rust_htslib::errors::Error::BamAuxTagNotFound);
        assert_eq!(record.aux(b"PS").err().unwrap(), rust_htslib::errors::Error::BamAuxTagNotFound);
    }

    #[test]
    fn test_pg_id_first_run_uses_base_id() {
        let version = "1.7.0-abcdef1";
        let base_id = format!("hiphase-v{version}");
        let header = append_and_roundtrip(&[], version);
        let ids = pg_ids_in_header(&header);
        assert_eq!(ids.len(), 1);
        assert!(ids.contains(&base_id));
    }

    #[test]
    fn test_pg_id_unrelated_prior_keeps_base_id() {
        let version = "1.7.0-abcdef1";
        let base_id = format!("hiphase-v{version}");
        let header = append_and_roundtrip(&["pbmm2", "hiphase-v1.6.0-oldhash"], version);
        let ids = pg_ids_in_header(&header);
        assert_eq!(ids.len(), 3);
        assert!(ids.contains(&"pbmm2".to_string()));
        assert!(ids.contains(&"hiphase-v1.6.0-oldhash".to_string()));
        assert!(ids.contains(&base_id));
        assert!(!ids.contains(&format!("{base_id}.1")));
    }

    #[test]
    fn test_pg_id_duplicate_gets_dot_one_suffix() {
        let version = "1.7.0-abcdef1";
        let base_id = format!("hiphase-v{version}");
        let header = append_and_roundtrip(&[&base_id], version);
        let ids = pg_ids_in_header(&header);
        assert_eq!(ids.len(), 2);
        assert!(ids.contains(&base_id));
        assert!(ids.contains(&format!("{base_id}.1")));
    }

    #[test]
    fn test_pg_id_third_run_gets_dot_two_suffix() {
        let version = "1.7.0-abcdef1";
        let base_id = format!("hiphase-v{version}");
        let prior_owned = [base_id.clone(), format!("{base_id}.1")];
        let prior_refs: Vec<&str> = prior_owned.iter().map(|s| s.as_str()).collect();
        let header = append_and_roundtrip(&prior_refs, version);
        let ids = pg_ids_in_header(&header);
        assert_eq!(ids.len(), 3);
        assert!(ids.contains(&base_id));
        assert!(ids.contains(&format!("{base_id}.1")));
        assert!(ids.contains(&format!("{base_id}.2")));
    }

    #[test]
    fn test_collect_pg_ids_ignores_co_lines() {
        let mut header = header_with_pg_ids(&["hiphase-v1.0"]);
        header.push_comment(b"ID:not-a-pg-id");
        let header_view = HeaderView::from_header(&header);
        let ids = collect_pg_ids(&header_view);
        assert_eq!(ids.len(), 1);
        assert!(ids.contains("hiphase-v1.0"));
    }
}