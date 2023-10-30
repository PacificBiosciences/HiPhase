
use crate::{block_gen::PhaseBlock, data_types::variants::VariantType};
use crate::data_types::read_segments::ReadSegment;
use crate::data_types::variants::Variant;
use crate::writers::phase_stats::PhaseStats;

use bio::data_structures::interval_tree::IntervalTree;
use log::{debug,trace};
use priority_queue::PriorityQueue;
use std::cmp::Reverse;

/// A node in the A* search tree.
#[derive(Eq,Hash,PartialEq)]
struct AstarNode {
    /// The node index
    node_index: u64,
    /// The cost that is fixed that no longer needs to be re-computed. Corresponds to reads that overlap early parts of the haplotypes.
    frozen_cost: u64,
    /// The cost that needs to be re-computed, corresponds to reads that partially overlap this solution, but have more variants.
    fluid_cost: u64,
    /// An estimate of the remaining cost to extend this node to full.
    heuristic_cost: u64,
    /// The first haplotype in this node's solution.
    h1: Vec<u8>,
    /// The second haplotype in this node's solution.  It can be identically to h1, but usually is not.
    h2: Vec<u8>,
    /// The number of heterozygous results in this node's solution.  I.e. sum(h1[x] != h2[x])
    num_hets: u64
}

impl std::fmt::Debug for AstarNode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("AstarNode")
            .field("frozen_cost", &self.frozen_cost)
            .field("fluid_cost", &self.fluid_cost)
            .field("heuristic_cost", &self.heuristic_cost)
            .field("hap.len()", &self.h1.len())
            .finish()
    }
}

impl AstarNode {
    /// Returns a new empty haplotype node with a heuristic cost.
    /// This should really only be used for a root node.
    /// # Arguments
    /// * `max_heuristic` - the estimate cost for the full phase block
    pub fn new(max_heuristic: u64) -> AstarNode {
        AstarNode {
            node_index: 0,
            frozen_cost: 0,
            fluid_cost: 0,
            heuristic_cost: max_heuristic,
            h1: Default::default(),
            h2: Default::default(),
            num_hets: 0
        }
    }

    /// This will create a newly extended node from haplotypes.
    /// Heuristic cost must be provided, but actual cost will be calculated from the reads.
    /// # Arguments
    /// * `node_index` - the index to use for this node, generally the index = the order of generating/encountering nodes
    /// * `parent_node` - the parent node used to spawn this one in our search tree space
    /// * `allele1` - the first allele, gets appended to haplotype 1 from parent
    /// * `allele2` - the second allele, gets appened to haplotype 2 from parent
    /// * `heuristic_cost` - the estimated cost of adding remaining variants to the haplotypes
    /// * `read_segments` - all the reads that we need to evaluate the _actual_ cost so far
    /// * `hap_offset` - the offset of the haplotype relative to read starts, only needed if solving subproblems
    pub fn new_extended_node(
        node_index: u64,
        parent_node: &AstarNode,
        allele1: u8, 
        allele2: u8,
        heuristic_cost: u64, 
        read_segments: &IntervalTree<usize, ReadSegment>,
        hap_offset: usize
    ) -> AstarNode {
        // make sure we didn't goof
        let mut h1: Vec<u8> = parent_node.get_h1().to_vec();
        h1.push(allele1);
        let mut h2: Vec<u8> = parent_node.get_h2().to_vec();
        h2.push(allele2);
        assert_eq!(h1.len(), h2.len());
        
        // num hets = parent hets + 0/1 depending on allelic extensions
        let num_hets = parent_node.get_num_hets()+(if allele1 == allele2 { 0 } else { 1 });

        // copy the frozen cost and initial fluid to 0, these will get added to below
        let mut frozen_cost = parent_node.get_frozen_cost();
        let mut fluid_cost = 0;
        let hap_len = h1.len()+hap_offset;
        for rs_interval in read_segments.find(hap_len-1..hap_len) {
            // calculate the cost of this segment with this phase block so far
            let rs = rs_interval.data();
            let rs_cost = std::cmp::min(
                rs.score_partial_haplotype(&h1[..], hap_offset),
                rs.score_partial_haplotype(&h2[..], hap_offset)
            );

            // determine where that cost gets assigned
            if rs.last_allele() < hap_len {
                // this one is frozen because the last allele was just added
                frozen_cost += rs_cost;
            } else {
                // there are more alleles for this segment, so it's liquid still
                fluid_cost += rs_cost;
            }
        }
        
        AstarNode {
            node_index,
            frozen_cost,
            fluid_cost,
            heuristic_cost,
            h1,
            h2,
            num_hets
        }
    }

    pub fn get_frozen_cost(&self) -> u64 {
        self.frozen_cost
    }

    /// Returns the combined frozen, fluid, and heuristic cost of this node.
    pub fn get_total_cost(&self) -> u64 {
        self.frozen_cost + self.fluid_cost + self.heuristic_cost
    }

    /// Priority is ranked by minimum total cost -> max number of hets -> earliest node index
    pub fn get_priority(&self) -> (Reverse<u64>, u64, Reverse<u64>) {
        (Reverse(self.get_total_cost()), self.num_hets, Reverse(self.node_index))
    }

    /// Returns a priority where cost is 0, this is primarily to trigger forced pruning
    pub fn get_cleared_priority(&self) -> (Reverse<u64>, u64, Reverse<u64>) {
        (Reverse(0), self.num_hets, Reverse(self.node_index))
    }

    pub fn get_h1(&self) -> &[u8] {
        &self.h1[..]
    }

    pub fn get_h2(&self) -> &[u8] {
        &self.h2[..]
    }

    pub fn get_allele_count(&self) -> usize {
        self.h1.len()
    }

    #[allow(dead_code)]
    pub fn get_node_index(&self) -> u64 {
        self.node_index
    }

    pub fn get_num_hets(&self) -> u64 {
        self.num_hets
    }

    /// Returns true if the internal haplotypes are identical.
    /// Usually only happens when near the root of the tree.
    pub fn is_identical_haplotypes(&self) -> bool {
        self.h1 == self.h2
    }
}

/// Struct for tracking the length of haplotypes in our queue that are greater than some threshold.
/// This allows us to track how long a queue actually is when we know we will ignore events smaller than the threshold.
/// Currently, it does not contain the queue itself, which may be worth doing in the long term.
struct PQueueHapTracker {
    /// The count of each haplotype size in the queue
    length_counts: Vec<usize>,
    /// The total number of haplotypes in the queue with length >= threshold
    total_count: usize,
    /// The minimum threshold for a haplotype to count
    threshold: usize
}

impl PQueueHapTracker {
    /// Creates a new tracker with a given maximum haplotype size.
    /// # Arguments
    /// * `max_hap_length` - the maximum size of the haplotypes tracked
    pub fn new(max_hap_length: usize) -> PQueueHapTracker {
        PQueueHapTracker {
            length_counts: vec![0; max_hap_length+1],
            total_count: 0,
            threshold: 0
        }
    }

    /// Adds a haplotype length to our tracker
    /// # Arguments
    /// * `value` - the length of the haplotype getting tracked
    pub fn add_hap(&mut self, value: usize) {
        self.length_counts[value] += 1;
        if value >= self.threshold {
            self.total_count += 1;
        }
    }

    /// Removes a haplotype length from the tracker
    /// # Arguments
    /// * `value` - the length of the haplotype getting removed from tracking
    pub fn remove_hap(&mut self, value: usize) {
        assert!(self.length_counts[value] > 0);
        self.length_counts[value] -= 1;
        if value >= self.threshold {
            assert!(self.total_count > 0);
            self.total_count -= 1;
        }
    }

    /// Increased the threshold of what is included in our total count
    /// # Arguments
    /// * `new_threshold` - the new minimum threshold to track, must be >= current threshold
    pub fn increase_threshold(&mut self, new_threshold: usize) {
        assert!(new_threshold >= self.threshold);
        trace!("increase_threshold => {}, size = {}", self.threshold, self.total_count);
        for t in self.threshold..new_threshold {
            self.total_count -= self.length_counts[t];
        }
        self.threshold = new_threshold;
        trace!("increase_threshold => {}, size = {}", self.threshold, self.total_count);
    }

    /// Returns the total number of haplotypes in the queue with length >= the internal threshold.
    pub fn len(&self) -> usize {
        self.total_count
    }
}

/// This calculates the heuristic estimates to be used by the full A* algorithm.
/// Given N alleles to phase, this will return a Vec with N+1 values, such that V[x] = H(x) where x is the number of set haplotype values.
/// This array is monotonically decreasing and always ends with 0.
/// Also returns a second boolean array indicating a variant is "problematic" and should be ignored (currently disabled, controlled via `bad_variants_enabled` constant).
/// This implementation uses the A* algorithm to generate the heuristic.
/// See `astar_subsolver(...)` for details.
/// # Arguments
/// * `num_variants` - the number of variants in the phase block; run-time grows linearly as the number of variants increases
/// * `max_segment_size` - the maximum number of variants to use when calculating sub-block heuristics; run-time grows a least linearly with this value
/// * `read_segments` - the reads to use when calculating the heuristics costs; run-time grows linearly with the length of this data
/// * `min_queue_size` - the minimum length of the queue
/// * `queue_increment` - the length that the queue grows as more variants are added to the solution
/// * `opt_bad_variants` - an optional set of "bad" or "ignored" variants that we should just ignore from the start
fn calculate_astar_heuristic(
    num_variants: usize, max_segment_size: usize, read_segments: &IntervalTree<usize, ReadSegment>,
    min_queue_size: usize, queue_increment: usize, opt_bad_variants: Option<Vec<bool>>
) -> (Vec<u64>, Vec<bool>) {
    assert!(max_segment_size >= 2);
    // an extra slot is included because it makes some checks go away later
    let mut heuristics: Vec<u64> = vec![0; num_variants+1];
    let mut bad_variants: Vec<bool> = match opt_bad_variants {
        Some(obv) => {
            assert_eq!(obv.len(), num_variants);
            obv
        },
        None => vec![false; num_variants]
    };
    let bad_variants_enabled = false;
    let mut max_clip_size: usize = 1;
    for v_index in (0..num_variants).rev() {
        debug!("solving subproblem {}..{}", v_index, v_index+max_clip_size);
        let (max_estimate, solve_size): (u64, usize) = astar_subsolver(
            v_index, max_clip_size, read_segments, &heuristics[..], &bad_variants[..],
            min_queue_size / 10, queue_increment
        );
        assert!(solve_size >= max_clip_size.min(2));
        debug!("  estimate => {} (solution distance => {}/{})", max_estimate, solve_size, max_clip_size);

        // monotonically decreasing so it's either
        // 1) the cost of the next entry OR
        // 2) this cost of solving this problem + the cost of solving everything after this problem
        if bad_variants_enabled && solve_size < max_clip_size {
            // we couldn't successfully form a partial haplotype
            bad_variants[v_index] = true;
        }

        if bad_variants[v_index] {
            // this one was a bad variant that will get ignored (e.g., 0-cost), just copy the heuristic
            heuristics[v_index] = heuristics[v_index+1];
        } else {
            // we successfully navigated estimation from this point
            assert!(max_estimate >= heuristics[v_index+1]);
            heuristics[v_index] = max_estimate;
        }

        max_clip_size = (solve_size+1).min(max_segment_size);
    }

    (heuristics, bad_variants)
}

/// An unpruned A* solver that will attempt to find the best path for a sub-problem.
/// It uses the heuristic estimates downstream to calculate this ones largest estimate.
/// For problem_size, `p`, and problem_offset, `o`, it will calculate max(`best_path(o..o+x)+H[o+x]`, for all `0 <= x <= p`).
/// This is basically the worst-case combination of two partial solutions within this region.
/// There is a very small queue used by this solver and if it reaches a maximum capacity, it will exit early without fully estimating.
/// However, as the solver goes, it will calculate the maximum heuristic encountered so far even if `x` does not get all the way to `p`.
/// While this estimator may exit early, it is not pruned, so you are guaranteed to find the best heuristic up to the point it exits.
/// In practice, most of them fully reach the problem size except in the problematic areas.
/// Returns tuple `(max_cost, farthest_estimate)` where max_cost is the highest heuristic cost found and farthest_estimate is how far the A* algorithm made it through the subproblem.
/// # Arguments
/// * `problem_offset` - the offset into the problem we are solving
/// * `problem_size` - the length of the sub-problem we are solving, usually constant except at the tail
/// * `read_segments` - the reads that are used to measure costs
/// * `heuristic_costs` - the heuristic costs so far, everything in `problem_offset+1..` is assumed to be already populated
/// * `bad_variants` - the bad variants encountered so far, if one is true, it will basically be ignored
/// * `min_queue_size` - the minimum length of the queue
/// * `queue_increment` - the length that the queue grows as more variants are added to the solution
fn astar_subsolver(
    problem_offset: usize, problem_size: usize, read_segments: &IntervalTree<usize, ReadSegment>, 
    heuristic_costs: &[u64], bad_variants: &[bool], min_queue_size: usize, queue_increment: usize
) -> (u64, usize) {
    // now, the core looping algorithm
    let mut pqueue: PriorityQueue<AstarNode, (Reverse<u64>, u64, Reverse<u64>)> = PriorityQueue::new();
    let mut next_node_index: u64 = 1;

    // this heuristic _should_ always be 0, because we're trying to calculate it
    assert_eq!(heuristic_costs[problem_offset], 0);
    // initialize with our neighbor heuristic (heuristic_costs is 1 longer than necessary, so no check needed here)
    let initial_estimate = heuristic_costs[problem_offset+1];
    
    // initialize this subsolver with the estimate from the _next_ node; it has to be >= that value
    let initial_node = AstarNode::new(initial_estimate);
    let initial_priority = initial_node.get_priority();
    pqueue.push(initial_node, initial_priority);

    let mut next_expected: usize = 0;
    let mut max_cost_so_far: u64 = 0;
    
    // we want a base level queue length, but it also needs to grow *slightly* with the length of the problem
    let max_visits: usize = min_queue_size + queue_increment * problem_size;
    let mut nodes_visited: usize = 0;

    //we loop as long as the next entry is shorter than the problem size
    while pqueue.peek().unwrap().0.get_allele_count() < problem_size && nodes_visited < max_visits {
        let (top_node, _top_priority) = pqueue.pop().unwrap();
        let allele_count: usize = top_node.get_allele_count();
        nodes_visited += 1;
        
        if allele_count == next_expected {
            // debug!("B#{} ({}/{}, {:?} {}/{}) => {:?}", phase_block.get_block_index(), next_expected, num_variants, top_priority, hap_tracker.len(), pqueue.len(), top_node);
            max_cost_so_far = max_cost_so_far.max(top_node.get_total_cost());
            next_expected += 1;
        }
        
        if bad_variants[problem_offset+allele_count] {
            // the next variant we want to add to the haplotype is a bad variant, so skip it
            let new_node = AstarNode::new_extended_node(
                next_node_index,
                &top_node, 2, 2, 
                heuristic_costs[problem_offset+allele_count+1], 
                read_segments,
                problem_offset
            );
            next_node_index += 1;

            // the new node should have identical total cost
            assert_eq!(top_node.get_total_cost(), new_node.get_total_cost());

            // add it to the queue
            let new_priority = new_node.get_priority();
            pqueue.push(new_node, new_priority);
        } else {
            // we didn't exit, so we need to add all expansions of this allele
            let hap_order = [(0, 1), (1, 0), (0, 0), (1, 1)];
            for &(h1, h2) in hap_order.iter() {
                // we don't want to add both 0-1 and 1-0 if the haplotypes before are identical; it doubles our work
                // so skip 1-0 if the haplotypes in the node are identical
                if !(h1 == 1 && h2 == 0 && top_node.is_identical_haplotypes()) {
                    // generate a new node and add to the queue
                    let new_node = AstarNode::new_extended_node(
                        next_node_index,
                        &top_node, h1, h2,
                        heuristic_costs[problem_offset+allele_count+1], 
                        read_segments,
                        problem_offset
                    );
                    next_node_index += 1;

                    // add it to the queue
                    let new_priority = new_node.get_priority();
                    pqueue.push(new_node, new_priority);
                }
            }
        }
    }

    if pqueue.peek().unwrap().0.get_allele_count() == problem_size {
        // loop terminated because we reached our problem size
        let (top_node, _top_priority) = pqueue.peek().unwrap();//pqueue.pop().unwrap();
        max_cost_so_far = max_cost_so_far.max(top_node.get_total_cost());
        next_expected += 1;
    } else {
        // we exited early, so whatever max we found is what we have
    }

    (max_cost_so_far, next_expected-1)
}

/// A result for a phasing algorithm, assumes diploid solution currently.
pub struct AstarResult {
    /// The first haplotype in the solution.
    pub haplotype_1: Vec<u8>,
    /// The second haplotype in the solution.
    pub haplotype_2: Vec<u8>,
    /// Optional statistics from the problem
    pub statistics: PhaseStats
}

/// Returns a phasing result by performing an A* tree search algorithm to calculate the best phase solution for a phase block.
/// This algorithm currently uses a fixed heuristic based on the distance to the end.
/// It also has a pruning strategy that fixes the priority queue size based on the farthest explored node so far.
/// # Arguments
/// * `phase_block` - the phase block summary information
/// * `variants` - the variants in the block, primarily passed-through to solution
/// * `read_segments` - interval tree of the reads that serve as data points for the phasing algorithm
/// * `min_queue_size` - the minimum length of the queue
/// * `queue_increment` - the length that the queue grows as more variants are added to the solution
pub fn astar_solver(
    phase_block: &PhaseBlock, variants: &[Variant], read_segments: &IntervalTree<usize, ReadSegment>,
    min_queue_size: usize, queue_increment: usize
) -> AstarResult {
    // we use this a lot
    let num_variants: usize = variants.len();

    // this is a sanity check for now that verifies that all read segments have all ignored variants set to 3 in that position
    // there is technically a cost here, but it is negligible for our sanity while debugging
    for rse in read_segments.find(0..usize::MAX) {
        let segment = rse.data();
        let alleles = segment.alleles();
        for (var_index, v) in variants.iter().enumerate() {
            if v.is_ignored() {
                assert!(alleles[var_index] == 3);
            }
        }
    }

    // add all the bad variants to this list that will seed the heuristic calculation list
    let bad_variants: Vec<bool> = variants.iter()
        .map(|v| v.is_ignored())
        .collect();
    let opt_bad_variants: Option<Vec<bool>> = Some(bad_variants);
    
    // our queue has a flat size + a buffer for each variant encountered so far
    let mut curr_queue_size_threshold: usize = min_queue_size;
    let full_prune_enabled: bool = true; // this will make sure the pqueue.len() field stays below the `max_queue_size` for conserving memory
    
    // set max queue size to either 2 * the highest functional queue size OR a large constant, whichever is greater
    // the large constant help prevent over-use of the full pruning algorithm if the base queue values are small
    let max_queue_size: usize = (2 * (min_queue_size + queue_increment * num_variants)).max(10000);
    let mut min_progress: usize = 0;
    let mut pqueue: PriorityQueue<AstarNode, (Reverse<u64>, u64, Reverse<u64>)> = PriorityQueue::new();
    let mut hap_tracker: PQueueHapTracker = PQueueHapTracker::new(num_variants);
    let mut next_expected = 0;
    
    // calculate a heuristic by looking ahead 40 variants
    // TODO: do we want to make this a parameter at some point?
    let max_segment_size: usize = 40;
    let (heuristic_costs, bad_variants): (Vec<u64>, Vec<bool>) = calculate_astar_heuristic(
        num_variants, max_segment_size, read_segments, min_queue_size, queue_increment, opt_bad_variants
    );
    debug!("Heuristics(<={}): {:?}", max_segment_size, heuristic_costs);
    debug!("Bad variants: {:?}", bad_variants);

    // sanity check, make sure any ignore variants are flagged in bad_variants now
    for (i, v) in variants.iter().enumerate() {
        if v.is_ignored() {
            assert!(bad_variants[i]);
        }
    }

    // statistics we want to gather
    let mut num_pruned: u64 = 0;
    let estimated_cost: u64 = heuristic_costs[0];

    // now, the core looping algorithm
    let initial_node = AstarNode::new(heuristic_costs[0]);
    let initial_priority = initial_node.get_priority();
    pqueue.push(initial_node, initial_priority);
    hap_tracker.add_hap(0);
    let mut next_node_index: u64 = 1;
    
    //we loop as long as the next entry is shorter than the expected number of variants
    while pqueue.peek().unwrap().0.get_allele_count() < num_variants {
        let (top_node, top_priority) = pqueue.pop().unwrap();
        let allele_count: usize = top_node.get_allele_count();
        hap_tracker.remove_hap(allele_count);
        trace!("popped: {:?} => {:?}, {:?}, {:?}", top_priority, top_node, top_node.get_h1(), top_node.get_h2());
        if allele_count == next_expected {
            debug!("B#{} ({}/{}, {:?} {}/{}) => {:?}", phase_block.get_block_index(), next_expected, num_variants, top_priority, hap_tracker.len(), pqueue.len(), top_node);
            next_expected += 1;
            if num_pruned == 0 {
                curr_queue_size_threshold += queue_increment;
                assert_eq!(curr_queue_size_threshold, min_queue_size + queue_increment * next_expected);
            }
        }

        // if a node has fewer alleles than our minimum progression, it gets cut
        if allele_count < min_progress {
            //println!("Pruning {:?} {:?}", top_priority, top_node);
            if num_pruned == 0 {
                // first time we do this, we need to clear our queue size back to the minimum
                curr_queue_size_threshold = min_queue_size;
            }
            num_pruned += 1;
            continue;
        }

        if bad_variants[allele_count] {
            // we are skipping this variant, generate a new node with unassigned values and add it back to the queue
            let new_node = AstarNode::new_extended_node(
                next_node_index,
                &top_node, 2, 2,
                heuristic_costs[allele_count+1], 
                read_segments,
                0
            );
            next_node_index += 1;

            let new_priority = new_node.get_priority();
            assert_eq!(top_node.get_total_cost(), new_node.get_total_cost());
            pqueue.push(new_node, new_priority);
            hap_tracker.add_hap(allele_count+1);
        } else {
            // we didn't exit, so we need to add all expansions of this allele
            // these are ordered such that heterozygous options come first
            let hap_order = [(0, 1), (1, 0), (0, 0), (1, 1)];
            for &(h1, h2) in hap_order.iter() {
                // we don't want to add both 0-1 and 1-0 if the haplotypes before are identical; it doubles our work
                // so skip 1-0 if the haplotypes in the node are identical
                if !(h1 == 1 && h2 == 0 && top_node.is_identical_haplotypes()) {
                    // generate a new node and add to the queue
                    let new_node = AstarNode::new_extended_node(
                        next_node_index,
                        &top_node, h1, h2,
                        heuristic_costs[allele_count+1], 
                        read_segments,
                        0
                    );
                    next_node_index += 1;

                    let new_priority = new_node.get_priority();
                    trace!("Pushing {:?}, {:?}, {:?}", new_node, new_node.get_h1(), new_node.get_h2());
                    pqueue.push(new_node, new_priority);
                    hap_tracker.add_hap(allele_count+1);
                }
            }
        }

        // check if we need to increase our minimum progression to prune off some nodes in the queue
        if hap_tracker.len() > curr_queue_size_threshold && min_progress < next_expected {
            min_progress += 1;
            debug!("B#{} min_progress={}", phase_block.get_block_index(), min_progress);
            hap_tracker.increase_threshold(min_progress);

            // check if our literal queue is holding too much data
            if full_prune_enabled && pqueue.len() > max_queue_size {
                // this means there are a lot of short haps we will eventually prune but haven't yet due to relative high cost
                // mark them as the "best" priority so they get deleted right away in the next few loops
                // this can be expensive, and has no functional benefit other than freeing memory
                // worst case HG001 test without this got up to 60GB
                let mut prune_count: usize = 0;
                for (node, priority) in pqueue.iter_mut() {
                    if node.get_allele_count() < min_progress {
                        *priority = node.get_cleared_priority();
                        prune_count += 1;
                    }
                }
                // when the iter_mut() is released, the pqueue automatically re-prioritizes itself, magic!
                debug!("B#{} Full prune triggered for {} nodes.", phase_block.get_block_index(), prune_count);
            }
        }
    }

    let (top_node, top_priority) = pqueue.pop().unwrap();
    let allele_count: usize = top_node.get_allele_count();
    hap_tracker.remove_hap(allele_count);
    if allele_count == num_variants {
        // successful full solution
        debug!("B#{} ({}/{}, {:?} {}/{}) => {:?}", phase_block.get_block_index(), next_expected, num_variants, top_priority, hap_tracker.len(), pqueue.len(), top_node);
        let haplotype_1 = top_node.get_h1().to_vec();
        let haplotype_2 = top_node.get_h2().to_vec();
        let actual_cost: u64 = top_node.get_total_cost();
        
        // gather stats on how many variants were phased, skipped, or homozygous in this block
        let mut phased_variants = 0;
        let mut phased_snvs = 0;
        let mut homozygous_variants = 0;
        let mut skipped_variants = 0;
        for (i, (&h1, &h2)) in haplotype_1.iter().zip(haplotype_2.iter()).enumerate() {
            if h1 != h2 {
                phased_variants += 1;
                if variants[i].get_type() == VariantType::Snv {
                    phased_snvs += 1;
                }
            } else if h1 == 2 {
                // they are both equal to 2
                skipped_variants += 1;
            } else {
                homozygous_variants += 1;
            }
        }
        debug!("B#{} phased: {}, homozygous: {}, skipped: {}", phase_block.get_block_index(), phased_variants, homozygous_variants, skipped_variants);

        let statistics: PhaseStats = PhaseStats::astar_new(
            num_pruned, estimated_cost, actual_cost,
            phased_variants, phased_snvs, homozygous_variants, skipped_variants
        );

        // send it all back
        AstarResult {
            haplotype_1,
            haplotype_2,
            statistics
        }
    } else {
        // given our current algorithm setup, this should never happen; we will always find _locally_ optimal paths through the tree
        panic!("B#{} failed to find solution; ({}/{}), {:?} {} => {:?}", phase_block.get_block_index(), next_expected, num_variants, top_priority, pqueue.len(), top_node);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    /// Returns two full-length reads: the first is all 0s with qual 2 and a second is all 1s with qual 3
    /// # Arguments
    /// * `num_alleles` - the length of each read
    fn get_simple_reads(num_alleles: usize) -> IntervalTree<usize, ReadSegment> {
        let rs1 = ReadSegment::new(
            "read_name".to_string(),
            vec![0; num_alleles],
            vec![2; num_alleles]
        );
        let rs2 = ReadSegment::new(
            "read_name_2".to_string(),
            vec![1; num_alleles],
            vec![3; num_alleles]
        );
        let seg_vec = vec![rs1, rs2];

        let mut read_segments: IntervalTree<usize, ReadSegment> = IntervalTree::new();
        for rs in seg_vec.into_iter() {
            read_segments.insert(rs.get_range(), rs);
        }
        read_segments
    }

    #[test]
    fn test_astarnode() {
        // all 0 read with quals = 2 and an all 1 read with quals = 3
        let num_alleles = 4;
        
        // it doesn't matter that these are wrong for testing
        let heuristic_costs: Vec<u64> = (0..(num_alleles+1)).map(|i| (num_alleles - i) as u64).collect();
        let hap_offset = 0;
        let read_segments = get_simple_reads(num_alleles);

        // test an all 0-hom mode
        let mut current_node = AstarNode::new(heuristic_costs[0]);
        for i in 0..num_alleles {
            let node_index: u64 = (i as u64)+1;
            let num_hets: u64 = 0;
            let next_node = AstarNode::new_extended_node(
                node_index,
                &current_node, 0, 0,
                heuristic_costs[i+1],
                &read_segments,
                hap_offset
            );

            let expected_cost = heuristic_costs[i+1] + 3 * node_index;
            assert_eq!(next_node.get_node_index(), node_index);
            
            let expected_frozen = if i == num_alleles - 1 {
                3 * num_alleles as u64 
            } else {
                0
            };
            assert_eq!(next_node.get_frozen_cost(), expected_frozen);
            assert_eq!(next_node.get_total_cost(), expected_cost);
            assert_eq!(next_node.get_priority(), (Reverse(expected_cost), num_hets, Reverse(node_index)));
            assert_eq!(next_node.get_cleared_priority(), (Reverse(0), num_hets, Reverse(node_index)));
            assert_eq!(next_node.get_h1().to_vec(), vec![0; i+1]);
            assert_eq!(next_node.get_h2().to_vec(), vec![0; i+1]);
            assert_eq!(next_node.get_allele_count(), i+1);
            assert_eq!(next_node.get_num_hets(), num_hets);
            assert!(next_node.is_identical_haplotypes());
            current_node = next_node;
        }

        // test an all het (0|1) mode
        let mut current_node = AstarNode::new(heuristic_costs[0]);
        for i in 0..num_alleles {
            let node_index: u64 = (i as u64)+1;
            let num_hets: u64 = node_index;
            let next_node = AstarNode::new_extended_node(
                node_index,
                &current_node, 0, 1,
                heuristic_costs[i+1],
                &read_segments,
                hap_offset
            );

            let expected_cost = heuristic_costs[i+1];
            assert_eq!(next_node.get_node_index(), node_index);
            
            let expected_frozen = 0;
            assert_eq!(next_node.get_frozen_cost(), expected_frozen);
            assert_eq!(next_node.get_total_cost(), expected_cost);
            assert_eq!(next_node.get_priority(), (Reverse(expected_cost), num_hets, Reverse(node_index)));
            assert_eq!(next_node.get_cleared_priority(), (Reverse(0), num_hets, Reverse(node_index)));
            assert_eq!(next_node.get_h1().to_vec(), vec![0; i+1]);
            assert_eq!(next_node.get_h2().to_vec(), vec![1; i+1]);
            assert_eq!(next_node.get_allele_count(), i+1);
            assert_eq!(next_node.get_num_hets(), num_hets);
            assert!(!next_node.is_identical_haplotypes());
            current_node = next_node;
        }

        // test an all 1-hom mode
        let mut current_node = AstarNode::new(heuristic_costs[0]);
        for i in 0..num_alleles {
            let node_index: u64 = (i as u64)+1;
            let num_hets: u64 = 0;
            let next_node = AstarNode::new_extended_node(
                node_index,
                &current_node, 1, 1,
                heuristic_costs[i+1],
                &read_segments,
                hap_offset
            );

            let expected_cost = heuristic_costs[i+1] + 2 * node_index;
            assert_eq!(next_node.get_node_index(), node_index);
            
            let expected_frozen = if i == num_alleles - 1 {
                2 * num_alleles as u64 
            } else {
                0
            };
            assert_eq!(next_node.get_frozen_cost(), expected_frozen);
            assert_eq!(next_node.get_total_cost(), expected_cost);
            assert_eq!(next_node.get_priority(), (Reverse(expected_cost), num_hets, Reverse(node_index)));
            assert_eq!(next_node.get_cleared_priority(), (Reverse(0), num_hets, Reverse(node_index)));
            assert_eq!(next_node.get_h1().to_vec(), vec![1; i+1]);
            assert_eq!(next_node.get_h2().to_vec(), vec![1; i+1]);
            assert_eq!(next_node.get_allele_count(), i+1);
            assert_eq!(next_node.get_num_hets(), num_hets);
            assert!(next_node.is_identical_haplotypes());
            current_node = next_node;
        }
    }

    #[test]
    fn test_pqueuehaptracker() {
        let mut tracker: PQueueHapTracker = PQueueHapTracker::new(10);
        for i in 0..11 {
            tracker.add_hap(i);
        }
        // make sure length matches everything so far
        assert_eq!(tracker.len(), 11);
        
        // try basic subtraction
        tracker.remove_hap(3);
        assert_eq!(tracker.len(), 10);

        // try threshold increasing
        tracker.increase_threshold(4);
        assert_eq!(tracker.len(), 7);
        
        // make sure removal of pre-threshold data doesn't change our len()
        for i in 0..3 {
            tracker.remove_hap(i);
            assert_eq!(tracker.len(), 7);
        }

        // make sure adding pre-threshold data doeesn't change our len()
        tracker.add_hap(0);
        assert_eq!(tracker.len(), 7);

        // make sure increasing threshold to the same value does nothing
        tracker.increase_threshold(4);
        assert_eq!(tracker.len(), 7);
    }
}