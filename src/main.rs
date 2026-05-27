
use hiphase::block_gen::{MultiPhaseBlockIterator, PhaseBlock, PhaseBlockIterator, get_sample_bams, get_vcf_samples};
use hiphase::cli::{Settings,check_settings,get_raw_settings};
use hiphase::data_types::reference_genome::ReferenceGenome;
use hiphase::data_types::variants::{VariantType, Zygosity};
use hiphase::phaser::{HaplotagResult, PhaseResult, solve_block, create_unphased_result};
use hiphase::writers::block_stats::BlockStatsCollector;
use hiphase::writers::haplotag_writer::HaplotagWriter;
use hiphase::writers::ordered_bam_writer::OrderedBamWriter;
use hiphase::writers::ordered_vcf_writer::OrderedVcfWriter;
use hiphase::writers::phase_stats::StatsWriter;
use hiphase::writers::vcf_util::build_bcf_index;

use log::{LevelFilter, debug, error, info, warn};
use rustc_hash::FxHashMap as HashMap;
use std::any::Any;
use std::path::PathBuf;
use std::sync::{Arc, mpsc};
use std::time::Instant;

/// Type alias for variant statistics from the block iterator
type VariantStats = HashMap<(String, String, VariantType, Zygosity), usize>;

/// Messages sent from the primary (block-generating) thread to the main thread
enum WorkerThreadMessage {
    /// A block count message, which will be sent after all blocks are submitted
    BlockCount(u64),
    /// A batch of phasing results to process
    BatchPhaseJobResult(Vec<(PhaseResult, HaplotagResult)>),
    /// Variant statistics collected from the block iterator
    VariantStats(Box<VariantStats>),
}

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

    // here's where the fun starts
    // generate blocks and sample-to-BAM mappings
    let (mut block_iterator, sample_to_bams, sample_to_output_bams) = 
        create_block_iterator(&cli_settings, &sample_names, &bam_thread_pool);

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

    // suppress repeated HTSlib warnings after the initial load
    suppress_htslib_warnings();

    //process the blocks (eventually in parallel)
    let start_time: Instant = Instant::now();
    let mut total_variants: u64 = 0;
    let mut results_received: u64 = 0;
    let mut variant_stats: Option<VariantStats> = None;
    
    // values related to printing
    const UPDATE_SPEED: u64 = 100;
    info!("Phase block generation starting...");

    if cli_settings.threads <= 1 {
        // single-threaded mode, so we just loop through the blocks sequentially and write results as we go
        let phasing_config = cli_settings.phasing_config();
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
                    &phasing_config
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

            if results_received.is_multiple_of(UPDATE_SPEED) {
                let time_so_far: f64 = start_time.elapsed().as_secs_f64();
                let blocks_per_sec: f64 = results_received as f64 / time_so_far;
                let variants_per_sec: f64 = total_variants as f64 / time_so_far;
                info!("Received results for {} phase blocks: {:.4} blocks/sec, {:.4} hets/sec, writer waiting on block {}", results_received, blocks_per_sec, variants_per_sec, vcf_writer.get_wait_block());
            }
        }

        // we need to get the variant stats from the block iterator because it will be dropped in the primary thread
        variant_stats = Some(block_iterator.variant_stats());
    } else {
        // we are going parallel, pull in rayon to handle the pool and threading
        // first, we need to drop the iterator because it will get recreated in the primary thread
        std::mem::drop(block_iterator);

        // clone the settings and sample names so we can pass them to the threads
        let cli_settings = cli_settings.clone();
        let sample_names = sample_names.clone();
        let arc_reference_genome = arc_reference_genome.clone();

        //set up job configuration
        info!("Starting job pool with {} threads...", cli_settings.threads);
        let mut jobs_queued: u64 = 0;
        
        // we need to set up the multiprocessing components now
        // a panic handler is provided so we can semi-cleanly exit the program if an unexpected panic occurs
        // if any of the threads use std::process::exit, the program will exit with the appropriate code still
        match rayon::ThreadPoolBuilder::new()
            .num_threads(cli_settings.threads as usize)
            .panic_handler(handle_thread_panic)
            .build_global() {
            Ok(()) => {},
            Err(e) => {
                error!("Error while building thread pool: {e}");
                std::process::exit(exitcode::OSERR);
            }
        };

        // channel for the primary thread to send messages back to the main thread
        let (tx, rx) = mpsc::channel::<WorkerThreadMessage>();
        let arc_phasing_config = Arc::new(cli_settings.phasing_config());
        let arc_cli_settings: Arc<Settings> = Arc::new(cli_settings.clone());
        let arc_sample_to_bams = Arc::new(sample_to_bams.clone());

        rayon::spawn(move || {
            // shared thread pool for bam IO
            let bam_thread_pool = match rust_htslib::tpool::ThreadPool::new(cli_settings.io_threads.unwrap() as u32) {
                Ok(btp) => btp,
                Err(e) => {
                    error!("Error while starting thread pool: {}", e);
                    std::process::exit(exitcode::IOERR);
                }
            };

            // generate blocks and sample-to-BAM mappings
            let (mut block_iterator, _sample_to_bams, _sample_to_output_bams) = 
                create_block_iterator(&arc_cli_settings, &sample_names, &bam_thread_pool);

            // we will collect the results here and send them back to the main thread in batches
            let batch_check_modulo = UPDATE_SPEED;
            let mut batch_problems: Vec<PhaseBlock> = vec![];

            // in an ideal world, we could use par_bridge here, but MultiPhaseBlockIterator does not implement Send (due to htslib)
            // if we ever decide to use something like noodles, we can revisit this, but it's not clear there's a benefit
            for (i, block_result) in block_iterator.by_ref()
                .enumerate().skip(skip_count).take(take_count) {
                let block = match block_result {
                    Ok(b) => b,
                    Err(e) => {
                        error!("Error while parsing VCF file: {}", e);
                        std::process::exit(exitcode::IOERR);
                    }
                };
                debug!("block {}: {:?} {}", i, block, block.bp_len());

                // update the job count and print an update if we're on the mod of our speed
                jobs_queued += 1;
                if jobs_queued.is_multiple_of(UPDATE_SPEED) {
                    info!("Generated {} phase blocks, latest block: {:?}", jobs_queued, block);
                }

                // check if we are phasing or short-circuiting
                if !block.unphased_block() && (phase_singletons || block.get_num_variants() > 1) {
                    // we are phasing, so clone all the Arcs and channels for the thread
                    let tx = tx.clone();
                    let arc_cli_settings = arc_cli_settings.clone();
                    let arc_phasing_config = arc_phasing_config.clone();
                    let arc_reference_genome = arc_reference_genome.clone();
                    let arc_sample_to_bams = arc_sample_to_bams.clone();

                    // spawn a thread to handle the phasing
                    rayon::spawn(move|| {
                        let sample_bams = arc_sample_to_bams.get(block.sample_name()).unwrap();

                        // dynamic errors cannot be sent via mpsc, so we need to handle errors here
                        let (phase_result, haplotag_result) = match solve_block(
                            &block,
                            &arc_cli_settings.vcf_filenames,
                            sample_bams,
                            &arc_reference_genome,
                            &arc_phasing_config
                        ) {
                            Ok(r) => r,
                            Err(e) => {
                                error!("Error while processing {:?}:", block);
                                error!("  {}", e);
                                std::process::exit(exitcode::SOFTWARE);
                            }
                        };

                        // send the result back to the main thread
                        tx.send(WorkerThreadMessage::BatchPhaseJobResult(
                            vec![(phase_result, haplotag_result)]
                        )).expect("channel will be there waiting for the pool");
                    });
                } else {
                    // short-circuiting, so we can just create the unphased result
                    // let (phase_result, haplotag_result) = create_unphased_result(&block);
                    batch_problems.push(block);
                }

                // check if we need to send a batch of simple problems to the thread
                // we do this in a separate thread because otherwise the main thread can get mutex-blocked on tx.send()
                if jobs_queued.is_multiple_of(batch_check_modulo) && !batch_problems.is_empty() {
                    debug!("Sending batch of {} problems to thread...", batch_problems.len());

                    // we need to send the batch of problems to the thread and clear create a new batch list
                    let thread_batches = batch_problems;
                    batch_problems = vec![];

                    // clone the channel for the thread
                    let tx = tx.clone();
                    rayon::spawn(move || {
                        // we need to create the unphased results for the batch
                        let batch_job_results: Vec<(PhaseResult, HaplotagResult)> = thread_batches.iter()
                            .map(|block| {
                                create_unphased_result(block)
                            })
                            .collect();

                        // send the batch of results back to the main thread
                        tx.send(WorkerThreadMessage::BatchPhaseJobResult(batch_job_results))
                            .expect("channel will be there waiting for the pool");
                    });
                }
            }

            // if there are any problems left, we need to send them to the thread
            if !batch_problems.is_empty() {
                debug!("Sending batch of {} problems to thread...", batch_problems.len());

                // we need to send the batch of problems to the thread and clear create a new batch list
                let thread_batches = batch_problems;

                // clone the channel for the thread
                let tx = tx.clone();

                rayon::spawn(move || {
                    // we need to create the unphased results for the batch
                    let batch_job_results: Vec<(PhaseResult, HaplotagResult)> = thread_batches.iter()
                        .map(|block| {
                            create_unphased_result(block)
                        })
                        .collect();

                    // send the batch of results back to the main thread
                    tx.send(WorkerThreadMessage::BatchPhaseJobResult(batch_job_results))
                        .expect("channel will be there waiting for the pool");
                });
            }
            
            // all jobs are submitted, send the block count and variant stats back to the main thread
            info!("Phase block generation complete, found {jobs_queued} phase blocks.");
            tx.send(WorkerThreadMessage::BlockCount(jobs_queued))
                .expect("channel will be there waiting for the pool");
            tx.send(WorkerThreadMessage::VariantStats(Box::new(block_iterator.variant_stats())))
                .expect("channel will be there waiting for the pool");
        });

        let mut jobs_queued: Option<u64> = None;
        while let Ok(message) = rx.recv() {
            // handle the message type we received
            match message {
                // a block count message, which will be sent after all blocks are submitted
                WorkerThreadMessage::BlockCount(count) => {
                    jobs_queued = Some(count);
                },
                // a batch of phase job results, which are primarily the short-circuiting results
                WorkerThreadMessage::BatchPhaseJobResult(batch_job_results) => {
                    // we need to process the batch of results
                    for (phase_result, haplotag_result) in batch_job_results.into_iter() {
                        // this is only for printing
                        total_variants += phase_result.phase_block.get_num_variants() as u64;
                        results_received += 1;

                        process_results(
                            phase_result, haplotag_result,
                            &mut stats_writer, &mut block_collector, &mut haplotag_writer,
                            &mut vcf_writer, &mut opt_bam_writers
                        );

                        // do an update if we're on the mod of our speed OR it's the last one for a thread
                        if results_received.is_multiple_of(UPDATE_SPEED) || (jobs_queued.unwrap_or(u64::MAX) - results_received) < cli_settings.threads as u64 {
                            let time_so_far: f64 = start_time.elapsed().as_secs_f64();
                            let blocks_per_sec: f64 = results_received as f64 / time_so_far;
                            let variants_per_sec: f64 = total_variants as f64 / time_so_far;
                            if let Some(jobs_queued) = jobs_queued {
                                info!("Received results for {} / {} phase blocks: {:.4} blocks/sec, {:.4} hets/sec, writer waiting on block {}", results_received, jobs_queued, blocks_per_sec, variants_per_sec, vcf_writer.get_wait_block());
                            } else {
                                info!("Received results for {} phase blocks: {:.4} blocks/sec, {:.4} hets/sec, writer waiting on block {}", results_received, blocks_per_sec, variants_per_sec, vcf_writer.get_wait_block());
                            }
                        }
                    }
                },
                // a variant stats message, which will be sent after all jobs are submitted and just goes into summary stats
                WorkerThreadMessage::VariantStats(stats) => {
                    // save the variant stats for later use
                    variant_stats = Some(*stats);
                }
            };
        }

        // the rx closed, so we should have all the results now
        if let Some(jobs_queued) = jobs_queued {
            if results_received != jobs_queued {
                error!("Results received count does not match jobs queued count, cannot finalize output files");
                std::process::exit(exitcode::SOFTWARE);
            }
        } else {
            error!("Jobs queued count was not received, cannot finalize output files");
            std::process::exit(exitcode::SOFTWARE);
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
    info!("Indexing output VCF files:");
    for vcf_fn in cli_settings.output_vcf_filenames.iter() {
        info!("\tIndexing {:?}...", vcf_fn);
        match build_bcf_index(vcf_fn, None, cli_settings.threads as u32, !cli_settings.csi_index) {
            Ok(()) => {},
            Err(e) => {
                error!("Error while building index for {:?}: {}", vcf_fn, e);
                std::process::exit(exitcode::IOERR);
            }
        };
    }
    info!("Finished indexing all output VCF files.");

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
        info!("Indexing output BAM files:");
        for bam_fn in cli_settings.output_bam_filenames.iter() {
            info!("\tIndexing {:?}...", bam_fn);
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
                Ok(()) => {},
                Err(e) => {
                    error!("Error while building index for {:?}: {}", bam_fn, e);
                    std::process::exit(exitcode::IOERR);
                }
            };
        }
        info!("Finished indexing all output BAM files.");
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
        let variant_stats = match variant_stats {
            Some(stats) => stats,
            None => {
                error!("Variant stats were not collected, cannot write summary statistics");
                std::process::exit(exitcode::IOERR);
            }
        };

        match block_collector.write_block_stats(
            &sample_names, filename, &arc_reference_genome, 
            variant_stats
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

/// A tuple mostly to make clippy happy
/// # Fields
/// * `block_iterator` - the block iterator for the combined samples
/// * `sample_to_bams` - a mapping of sample names to the BAM files used for that sample
/// * `sample_to_output_bams` - a mapping of sample names to the BAM files to write to
type BlockIteratorFields = (MultiPhaseBlockIterator, HashMap<String, Vec<PathBuf>>, HashMap<String, Vec<PathBuf>>);

/// Creates a MultiPhaseBlockIterator along with sample-to-BAM mappings.
/// This function *must* be deterministic, otherwise there could be a desync in the block iterator and the writers.
/// # Arguments
/// * `cli_settings` - the CLI settings
/// * `sample_names` - the sample names to process
/// * `bam_thread_pool` - the thread pool for BAM IO
fn create_block_iterator(
    cli_settings: &Settings,
    sample_names: &[String],
    bam_thread_pool: &rust_htslib::tpool::ThreadPool,
) -> BlockIteratorFields {
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
            bam_thread_pool
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
    let block_iterator: MultiPhaseBlockIterator = match MultiPhaseBlockIterator::new(block_iterators) {
        Ok(mpbi) => mpbi,
        Err(e) => {
            error!("Error during phase block iterator creation: {}", e);
            std::process::exit(exitcode::IOERR);
        }
    };

    (block_iterator, sample_to_bams, sample_to_output_bams)
}

/// Suppress HTSlib warning-level log output for the remainder of the process.
/// Affects all threads (global `hts_verbose`); hides all HTSlib warnings, not only stale-index messages.
fn suppress_htslib_warnings() {
    unsafe {
        rust_htslib::htslib::hts_set_log_level(
            rust_htslib::htslib::htsLogLevel_HTS_LOG_ERROR,
        );
    }
}

/// Panic handler for the rayon thread pool that logs an error and exits the program
/// # Arguments
/// * `panic_info` - the panic information from the thread pool
fn handle_thread_panic(panic_info: Box<dyn Any + Send>) {
    // Try to extract a meaningful message from the panic payload
    let msg = if let Some(s) = panic_info.downcast_ref::<&str>() {
        s.to_string()
    } else if let Some(s) = panic_info.downcast_ref::<String>() {
        s.clone()
    } else {
        "Unknown panic".to_string()
    };

    error!("Panic detected in thread pool: {msg}");
    error!("Exiting program due to panic in thread pool");
    std::process::exit(exitcode::SOFTWARE);
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