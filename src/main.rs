
use hiphase::block_gen::{MultiPhaseBlockIterator, PhaseBlockIterator, get_vcf_samples, get_sample_bams};
use hiphase::cli::{Settings,check_settings,get_raw_settings};
use hiphase::data_types::reference_genome::ReferenceGenome;
use hiphase::phaser::{HaplotagResult, PhaseResult, solve_block, create_unphased_result};
use hiphase::writers::block_stats::BlockStatsCollector;
use hiphase::writers::haplotag_writer::HaplotagWriter;
use hiphase::writers::ordered_bam_writer::OrderedBamWriter;
use hiphase::writers::ordered_vcf_writer::OrderedVcfWriter;
use hiphase::writers::phase_stats::StatsWriter;
use hiphase::writers::vcf_util::build_bcf_index;

use log::{LevelFilter, debug, error, info, warn};
use rustc_hash::FxHashMap as HashMap;
use std::path::PathBuf;
use std::sync::{Arc, mpsc};
use std::time::Instant;
use threadpool::ThreadPool;

fn main() {
    // get the settings
    let settings: Settings = get_raw_settings();
    let filter_level: LevelFilter = match settings.verbosity {
        0 => LevelFilter::Info,
        1 => LevelFilter::Debug,
        _ => LevelFilter::Trace
    };

    // immediately setup logging first
    env_logger::builder()
        .format_timestamp_millis()
        .filter_level(filter_level)
        .init();
    
    // okay, now we can check all the other settings
    let cli_settings: Settings = check_settings(settings);

    // first we need to figure out which samples are getting phased
    let mut sample_names: Vec<String> = cli_settings.sample_names.clone();
    if sample_names.is_empty() {
        // no samples were provided, so add the first one encountered
        // we need to just infer that we're phasing the first one only for now
        let all_sample_names = match get_vcf_samples(&cli_settings.vcf_filenames[0]) {
            Ok(v) => v,
            Err(e) => {
                error!("Error during VCF sample name parsing: {}", e);
                std::process::exit(exitcode::IOERR);
            }
        };

        // some warnings as needed
        if all_sample_names.len() > 1 {
            warn!("Multi-sample VCF detected, but sample name was not provided.  Assuming name is {:?}.", all_sample_names[0]);
        } else {
            debug!("Single-sample VCF detected, but sample name was not provided.  Assuming name is {:?}.", all_sample_names[0]);
        }
        sample_names.push(all_sample_names[0].clone());
    }

    // if we are ignoring read groups, we need to verify only one sample is in use
    if cli_settings.ignore_read_groups && sample_names.len() > 1 {
        error!("Flag --ignore-read-groups cannot be used in conjuction with multiple sample names, either add read groups or run one sample name at a time.");
        std::process::exit(exitcode::USAGE);
    }

    // shared thread pool for bam IO
    let bam_thread_pool = match rust_htslib::tpool::ThreadPool::new(cli_settings.io_threads.unwrap() as u32) {
        Ok(btp) => btp,
        Err(e) => {
            error!("Error while starting thread pool: {}", e);
            std::process::exit(exitcode::IOERR);
        }
    };

    //here's where the fun starts
    //generate blocks
    let mut block_iterators: Vec<PhaseBlockIterator> = vec![];
    let mut all_used_bams = vec![];
    let mut sample_to_bams: HashMap<String, Vec<PathBuf>> = Default::default();
    let mut sample_to_output_bams: HashMap<String, Vec<PathBuf>> = Default::default();
    for sample_name in sample_names.iter() {
        // figure out which BAMS go with the given sample
        let (mut sample_bams, bam_indices) = if cli_settings.ignore_read_groups {
            // if we are ignoring read groups, then we use all bams (and all indices)
            (
                cli_settings.bam_filenames.clone(),
                (0..cli_settings.bam_filenames.len()).collect()
            )
        } else {
            match get_sample_bams(&cli_settings.bam_filenames, sample_name, &cli_settings.reference_filename) {
                Ok(sb) => sb,
                Err(e) => {
                    error!("Error during BAM read group parsing: {}", e);
                    std::process::exit(exitcode::IOERR);
                }
            }
        };
        sample_to_bams.insert(sample_name.clone(), sample_bams.clone());

        // make a phase block iterator using just the sample-specific bams
        let block_iterator: PhaseBlockIterator = match PhaseBlockIterator::new(
            &cli_settings.vcf_filenames,
            &sample_bams,
            &cli_settings.reference_filename,
            sample_name.clone(),
            cli_settings.min_variant_quality,
            cli_settings.min_mapping_quality,
            cli_settings.min_spanning_reads,
            !cli_settings.disable_supplemental_joins,
            &bam_thread_pool
        ) {
            Ok(bi) => bi,
            Err(e) => {
                error!("Error during file loading: {}", e);
                std::process::exit(exitcode::IOERR);
            }
        };

        // add the iterator to our list to put together
        block_iterators.push(block_iterator);

        // also save the used bams, we will check these soon
        all_used_bams.append(&mut sample_bams);

        // check if we need to save names for BAM writing
        if !cli_settings.output_bam_filenames.is_empty() {
            let mut sample_output_bams = vec![];
            for &b_index in bam_indices.iter() {
                sample_output_bams.push(cli_settings.output_bam_filenames[b_index].clone());
            }
            sample_to_output_bams.insert(sample_name.clone(), sample_output_bams);
        }
    }

    if cli_settings.bam_filenames.len() != all_used_bams.len() {
        let num_provided = cli_settings.bam_filenames.len();
        let num_used = all_used_bams.len();
        error!("User provided {} BAM files, but only {} matched samples for phasing", num_provided, num_used);
        error!("Please remove extra BAM files or add additional samples, BAMs matching phasing: {:?}", all_used_bams);
        std::process::exit(exitcode::IOERR);
    }

    // create our joint iterator
    let mut block_iterator: MultiPhaseBlockIterator = match MultiPhaseBlockIterator::new(block_iterators) {
        Ok(mpbi) => mpbi,
        Err(e) => {
            error!("Error during phase block iterator creation: {}", e);
            std::process::exit(exitcode::IOERR);
        }
    };

    // this writer will write "in-order" provided we correctly pass the ordering of data to it
    let mut vcf_writer: OrderedVcfWriter = match OrderedVcfWriter::new(
        &cli_settings.vcf_filenames,
        &cli_settings.output_vcf_filenames,
        cli_settings.min_variant_quality,
        &sample_names
    ) {
        Ok(vw) => vw,
        Err(e) => {
            error!("Error during VCF writer creation: {}", e);
            std::process::exit(exitcode::IOERR);
        }
    };

    // this write will write reads "in-order" provided we correctly pass the ordering of data to it
    let mut opt_bam_writers: Option<HashMap<String, OrderedBamWriter>> = if cli_settings.output_bam_filenames.is_empty() {
        None
    } else {
        let mut writer_map: HashMap<String, OrderedBamWriter> = Default::default();
        for sample_name in sample_names.iter() {
            let sample_bams = sample_to_bams.get(sample_name).unwrap();
            let sample_output_bams = sample_to_output_bams.get(sample_name).unwrap();
            writer_map.insert(
                sample_name.clone(), 
                match OrderedBamWriter::new(
                    sample_name.clone(),
                    &cli_settings.reference_filename,
                    sample_bams,
                    sample_output_bams,
                    &bam_thread_pool
                ) {
                    Ok(bw) => bw,
                    Err(e) => {
                        error!("Error during BAM writer creation: {}", e);
                        std::process::exit(exitcode::IOERR);
                    }
                }
            );
        }
        Some(writer_map)
    };

    // create our stats file also
    let mut stats_writer: Option<StatsWriter> = match cli_settings.stats_filename {
        Some(ref filename) => {
            match StatsWriter::new(filename) {
                Ok(sw) => Some(sw),
                Err(e) => {
                    error!("Error during statistics writer creation: {}", e);
                    std::process::exit(exitcode::IOERR);
                }
            }
        },
        None => None
    };

    // create our block stats collector
    let mut block_collector: BlockStatsCollector = BlockStatsCollector::new();

    let skip_count = cli_settings.skip_blocks;
    let take_count = cli_settings.take_blocks;
    let debug_run: bool = if skip_count != 0 || take_count != usize::MAX {
        warn!("Debug run detected, disabling file finalizing steps.");
        warn!("Blocks to skip: {}", skip_count);
        warn!("Blocks to process: {}", take_count);
        true
    } else {
        false
    };

    // create our haplotag file if necessary
    let mut haplotag_writer: Option<HaplotagWriter> = match cli_settings.haplotag_filename {
        Some(ref filename) => {
            match HaplotagWriter::new(filename) {
                Ok(hw) => Some(hw),
                Err(e) => {
                    error!("Error during haplotag writer creations: {}", e);
                    std::process::exit(exitcode::IOERR);
                }
            }
        },
        None => None
    };

    // controls whether singletons are deeply run, including haplotagging
    let phase_singletons: bool = cli_settings.phase_singletons;

    // get our reference genome if we have one
    let reference_genome: ReferenceGenome = match ReferenceGenome::from_fasta(&cli_settings.reference_filename) {
        Ok(rg) => rg,
        Err(e) => {
            error!("Error during reference loading: {}", e);
            std::process::exit(exitcode::IOERR);
        }
    };

    // check if the proper indexing is enabled; tbi/bai go up to 2**29 - 1 before bailing
    let max_chrom_len = reference_genome.contig_keys().iter()
        .map(|k| reference_genome.get_full_chromosome(k).len())
        .max()
        .unwrap_or_default();
    let csi_required = max_chrom_len > (2_usize.pow(29) - 1);
    if csi_required && !cli_settings.csi_index {
        error!("Output files will require .csi indexing ({max_chrom_len} > 2^29 - 1); use --csi-index to enable");
        std::process::exit(exitcode::USAGE);
    }

    // we have to do this because we need access to the reference genome later also
    let arc_reference_genome: Arc<ReferenceGenome> = Arc::new(reference_genome);

    //process the blocks (eventually in parallel)
    let start_time: Instant = Instant::now();
    let mut total_variants: u64 = 0;
    let mut results_received: u64 = 0;
    
    // values related to printing
    const UPDATE_SPEED: u64 = 100;
    info!("Phase block generation starting...");

    if cli_settings.threads <= 1 {
        for (i, block_result) in block_iterator.by_ref().enumerate().skip(skip_count).take(take_count) {
            let block = match block_result {
                Ok(b) => b,
                Err(e) => {
                    error!("Error while parsing VCF file: {}", e);
                    std::process::exit(exitcode::IOERR);
                }
            };
            debug!("block {}: {:?} {}", i, block, block.bp_len());

            // we likely need to separate out the phase result from the haplotag result
            let sample_bams = sample_to_bams.get(block.sample_name()).unwrap();
            let (phase_result, haplotag_result): (PhaseResult, HaplotagResult) = if !block.unphased_block() && (phase_singletons || block.get_num_variants() > 1) {
                match solve_block(
                    &block,
                    &cli_settings.vcf_filenames,
                    sample_bams,
                    &arc_reference_genome,
                    cli_settings.reference_buffer,
                    cli_settings.min_matched_alleles,
                    cli_settings.min_mapping_quality,
                    cli_settings.global_realign_cputime,
                    cli_settings.phase_min_queue_size,
                    cli_settings.phase_queue_increment,
                    cli_settings.wfa_prune_distance
                ) {
                    Ok(r) => r,
                    Err(e) => {
                        error!("Error while processing {:?}:", block);
                        error!("  {}", e);
                        std::process::exit(exitcode::SOFTWARE);
                    }
                }
            } else {
                create_unphased_result(&block)
            };

            // this is only for printing
            total_variants += phase_result.phase_block.get_num_variants() as u64;
            results_received += 1;

            process_results(
                phase_result, haplotag_result, 
                &mut stats_writer, &mut block_collector, &mut haplotag_writer,
                &mut vcf_writer, &mut opt_bam_writers,
            );

            if results_received % UPDATE_SPEED == 0 {
                let time_so_far: f64 = start_time.elapsed().as_secs_f64();
                let blocks_per_sec: f64 = results_received as f64 / time_so_far;
                let variants_per_sec: f64 = total_variants as f64 / time_so_far;
                info!("Received results for {} phase blocks: {:.4} blocks/sec, {:.4} hets/sec, writer waiting on block {}", results_received, blocks_per_sec, variants_per_sec, vcf_writer.get_wait_block());
            }
        }
    } else {
        //set up job configuration
        info!("Starting job pool with {} threads...", cli_settings.threads);
        let job_slots: u64 = 40 * cli_settings.threads as u64;
        let mut jobs_queued: u64 = 0;
        
        //we need to set up the multiprocessing components now
        let pool = ThreadPool::new(cli_settings.threads);
        let (tx, rx) = mpsc::channel();
        let arc_cli_settings: Arc<Settings> = Arc::new(cli_settings.clone());
        let arc_sample_to_bams = Arc::new(sample_to_bams.clone());

        for (i, block_result) in block_iterator.by_ref().enumerate().skip(skip_count).take(take_count) {
            // make sure no panics encountered so far
            if pool.panic_count() > 0 {
                error!("Panic detected in ThreadPool, check above for details.");
                std::process::exit(exitcode::SOFTWARE);
            }

            if jobs_queued - results_received >= job_slots {
                let (phase_result, haplotag_result): (PhaseResult, HaplotagResult) = rx.recv().unwrap();
                
                // this is only for printing
                total_variants += phase_result.phase_block.get_num_variants() as u64;
                results_received += 1;
                                
                process_results(
                    phase_result, haplotag_result, 
                    &mut stats_writer, &mut block_collector, &mut haplotag_writer,
                    &mut vcf_writer, &mut opt_bam_writers
                );

                if results_received % UPDATE_SPEED == 0 {
                    let time_so_far: f64 = start_time.elapsed().as_secs_f64();
                    let blocks_per_sec: f64 = results_received as f64 / time_so_far;
                    let variants_per_sec: f64 = total_variants as f64 / time_so_far;
                    info!("Received results for {} phase blocks: {:.4} blocks/sec, {:.4} hets/sec, writer waiting on block {}", results_received, blocks_per_sec, variants_per_sec, vcf_writer.get_wait_block());
                }
            }
            
            let block = match block_result {
                Ok(b) => b,
                Err(e) => {
                    error!("Error while parsing VCF file: {}", e);
                    std::process::exit(exitcode::IOERR);
                }
            };
            debug!("block {}: {:?} {}", i, block, block.bp_len());

            jobs_queued += 1;
            if jobs_queued % UPDATE_SPEED == 0 {
                info!("Generated {} phase blocks, latest block: {:?}", jobs_queued, block); 
            }

            if !block.unphased_block() && (phase_singletons || block.get_num_variants() > 1) {
                let tx = tx.clone();
                let arc_cli_settings = arc_cli_settings.clone();
                let arc_reference_genome = arc_reference_genome.clone();
                let arc_sample_to_bams = arc_sample_to_bams.clone();
                
                pool.execute(move|| {
                    let sample_bams = arc_sample_to_bams.get(block.sample_name()).unwrap();
                    // dynamic errors cannot be sent via mpsc, so we need to handle errors here
                    let all_results = match solve_block(
                        &block,
                        &arc_cli_settings.vcf_filenames,
                        sample_bams,
                        &arc_reference_genome,
                        arc_cli_settings.reference_buffer,
                        arc_cli_settings.min_matched_alleles,
                        arc_cli_settings.min_mapping_quality,
                        arc_cli_settings.global_realign_cputime,
                        arc_cli_settings.phase_min_queue_size,
                        arc_cli_settings.phase_queue_increment,
                        arc_cli_settings.wfa_prune_distance
                    ) {
                        Ok(r) => r,
                        Err(e) => {
                            error!("Error while processing {:?}:", block);
                            error!("  {}", e);
                            std::process::exit(exitcode::SOFTWARE);
                        }
                    };
                    tx.send(all_results).expect("channel will be there waiting for the pool");
                });
            } else {
                // this is a unphased block we can short-circuit here
                let (phase_result, haplotag_result): (PhaseResult, HaplotagResult) = 
                    create_unphased_result(&block);
                
                // this is only for printing
                total_variants += phase_result.phase_block.get_num_variants() as u64;
                results_received += 1;
                                
                process_results(
                    phase_result, haplotag_result, 
                    &mut stats_writer, &mut block_collector, &mut haplotag_writer,
                    &mut vcf_writer, &mut opt_bam_writers
                );

                if results_received % UPDATE_SPEED == 0 {
                    let time_so_far: f64 = start_time.elapsed().as_secs_f64();
                    let blocks_per_sec: f64 = results_received as f64 / time_so_far;
                    let variants_per_sec: f64 = total_variants as f64 / time_so_far;
                    info!("Received results for {} phase blocks: {:.4} blocks/sec, {:.4} hets/sec, writer waiting on block {}", results_received, blocks_per_sec, variants_per_sec, vcf_writer.get_wait_block());
                }
            }
        }

        while results_received < jobs_queued {
            // make sure no panics encountered so far
            // TODO: if we hit deadlocks from panics, we may need to add this with some sort of suitable timeout:
            //       https://doc.rust-lang.org/std/sync/mpsc/struct.Receiver.html#method.recv_timeout
            if pool.panic_count() > 0 {
                error!("Panic detected in ThreadPool, check above for details.");
                std::process::exit(exitcode::SOFTWARE);
            }

            let (phase_result, haplotag_result): (PhaseResult, HaplotagResult) = rx.recv().unwrap();
            
            // this is only for printing
            total_variants += phase_result.phase_block.get_num_variants() as u64;
            results_received += 1;
            
            process_results(
                phase_result, haplotag_result, 
                &mut stats_writer, &mut block_collector, &mut haplotag_writer,
                &mut vcf_writer, &mut opt_bam_writers
            );

            // do an update if we're on the mod of our speed OR it's the last one for a thread
            if results_received % UPDATE_SPEED == 0 || (jobs_queued - results_received) < cli_settings.threads as u64 {
                let time_so_far: f64 = start_time.elapsed().as_secs_f64();
                let blocks_per_sec: f64 = results_received as f64 / time_so_far;
                let variants_per_sec: f64 = total_variants as f64 / time_so_far;
                info!("Received results for {} / {} phase blocks: {:.4} blocks/sec, {:.4} hets/sec, writer waiting on block {}", results_received, jobs_queued, blocks_per_sec, variants_per_sec, vcf_writer.get_wait_block());
            }
        }
    }

    info!("All phase blocks analyzed, finalizing output files...");

    // if we are only doing partial files, this will not behave, so skip it
    if !debug_run {
        // we call this once at the end
        match vcf_writer.write_to_end_position() {
            Ok(()) => {},
            Err(e) => {
                error!("Error while finalizing VCF chromosomes: {}", e);
                std::process::exit(exitcode::IOERR);
            }
        };
    }

    // now we drop the VCF writer, this is to close out all the VCF files before indexing
    std::mem::drop(vcf_writer);
    for vcf_fn in cli_settings.output_vcf_filenames.iter() {
        match build_bcf_index(vcf_fn, None, cli_settings.threads as u32, !cli_settings.csi_index) {
            Ok(()) => {
                info!("Finished building index for {:?}.", vcf_fn);
            },
            Err(e) => {
                error!("Error while building index for {:?}: {}", vcf_fn, e);
                std::process::exit(exitcode::IOERR);
            }
        };
    }

    if let Some(bam_writers) = opt_bam_writers.as_mut() {
        // if we are only doing partial files, this will not behave, so skip it
        if !debug_run {
            for bam_writer in bam_writers.values_mut() {
                // first finalize whichever chromosome we were on
                match bam_writer.finalize_chromosome() {
                    Ok(()) => {},
                    Err(e) => {
                        error!("Error while finalizing BAM chromosomes: {}", e);
                        std::process::exit(exitcode::IOERR);
                    }
                };

                // copy reads from all the remaining chromosomes
                match bam_writer.copy_remaining_chromosomes() {
                    Ok(()) => {},
                    Err(e) => {
                        error!("Error while copying all remaining chromosomes: {}", e);
                        std::process::exit(exitcode::IOERR);
                    }
                };
            }
        }
        
        // now we need to drop the bam writer, this is to close out all the BAM files before indexing
        std::mem::drop(opt_bam_writers);

        // index the BAM files with .bai files
        for bam_fn in cli_settings.output_bam_filenames.iter() {
            let idx_type = if cli_settings.csi_index {
                rust_htslib::bam::index::Type::Csi(14)  
            } else {
                rust_htslib::bam::index::Type::Bai
            };
            match rust_htslib::bam::index::build(
                bam_fn,
                None,
                idx_type,
                cli_settings.threads as u32
            ) {
                Ok(()) => {
                    info!("Finished building index for {:?}.", bam_fn);
                },
                Err(e) => {
                    error!("Error while building index for {:?}: {}", bam_fn, e);
                    std::process::exit(exitcode::IOERR);
                }
            };
        }
    }

    if let Some(ref filename) = cli_settings.blocks_filename {
        // this will save all block information to a csv/tsv file
        info!("Saving all blocks to {:?}...", filename);
        match block_collector.write_blocks(filename) {
            Ok(()) => {},
            Err(e) => {
                error!("Error while writing blocks file: {}", e);
                std::process::exit(exitcode::IOERR);
            }
        };
    }

    if let Some(ref filename) = cli_settings.summary_filename {
        // this will save chromosome level stats to a csv/tsv file
        info!("Saving summary block statistics to {:?}...", filename);
        match block_collector.write_block_stats(
            &sample_names, filename, &arc_reference_genome, 
            block_iterator.variant_stats()
        ) {
            Ok(()) => {},
            Err(e) => {
                error!("Error while writing summary statistics file: {}", e);
                std::process::exit(exitcode::IOERR);
            }
        }
    }

    info!("All phase blocks finished successfully after {} seconds.", start_time.elapsed().as_secs_f64());
}

/// Sub-routine to make sure we are always consistently processing results in an identical manner
/// This is mostly because I got tired of forgetting to change things in 3 places
/// # Argument
/// * `phase_result` - the phasing result from our algorithm
/// * `haplotag_result` - the haplotag result, does nothing if we are not haplotagging
/// * `opt_stats_writer` - mutable, optional reference to our algorithm stats writer
/// * `block_collector` - mutable reference to the block stats collector
/// * `opt_haplotag_writer` - mutable, optional reference to our haplotag writer
/// * `vcf_writer` - mutable reference to our VCF writer
/// * `opt_bam_writers` - mutable, optional reference to the BAM writers for haplotagging
fn process_results(
    phase_result: PhaseResult, haplotag_result: HaplotagResult,
    opt_stats_writer: &mut Option<StatsWriter>, block_collector: &mut BlockStatsCollector, 
    opt_haplotag_writer: &mut Option<HaplotagWriter>,
    vcf_writer: &mut OrderedVcfWriter, opt_bam_writers: &mut Option<HashMap<String, OrderedBamWriter>>
) {
    // common debug statements
    debug!("block {} haplotypes:", phase_result.phase_block.get_block_index());
    debug!("{:?}", phase_result.haplotype_1);
    debug!("{:?}", phase_result.haplotype_2);

    // write the stats if we have both a writer and a stats block
    if let Some(stats_writer) = opt_stats_writer.as_mut() {
        match stats_writer.write_stats(&phase_result) {
            Ok(()) => {},
            Err(e) => {
                error!("Error while writing statistics file: {}", e);
                std::process::exit(exitcode::IOERR);
            }
        }
    };

    // save all the blocks here
    for sub_block in phase_result.sub_phase_blocks.iter() {
        block_collector.add_block(sub_block.clone());
    }
    block_collector.add_result(&phase_result);

    match vcf_writer.write_phase_block(phase_result) {
        Ok(()) => {},
        Err(e) => {
            error!("Error while saving phase block: {}", e);
            std::process::exit(exitcode::IOERR);
        }
    };

    if let Some(haplotag_writer) = opt_haplotag_writer.as_mut() {
        match haplotag_writer.write_block(&haplotag_result) {
            Ok(()) => {},
            Err(e) => {
                error!("Error while writing haplotag file: {}", e);
                std::process::exit(exitcode::IOERR);
            }
        };
    }

    if let Some(bam_writers) = opt_bam_writers.as_mut() {
        let sample_name = haplotag_result.phase_block.sample_name().to_string();
        let block_index = haplotag_result.phase_block.get_block_index();

        // send this block to the correct writer
        let bam_writer = bam_writers.get_mut(&sample_name).unwrap();
        match bam_writer.write_phase_block(haplotag_result) {
            Ok(()) => {},
            Err(e) => {
                error!("Error while saving haplotags: {}", e);
                std::process::exit(exitcode::IOERR);
            }
        };

        // now send the skip signal to the rest
        for (sn, bam_writer) in bam_writers.iter_mut() {
            if sn != &sample_name {
                match bam_writer.write_dummy_block(block_index) {
                    Ok(()) => {},
                    Err(e) => {
                        error!("Error while saving haplotags: {}", e);
                        std::process::exit(exitcode::IOERR);
                    }
                };
            }
        }
    }
}