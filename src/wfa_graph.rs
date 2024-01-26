
use crate::data_types::variants::Variant;

use bit_vec::BitVec;
#[allow(unused_imports)]
use log::{debug, trace, warn};
use priority_queue::PriorityQueue;
use simple_error::bail;
use std::cmp::Reverse;
use rustc_hash::FxHashMap as HashMap;

pub type NodeAlleleMap = HashMap<usize, Vec<(usize, u8)>>;

/// Contains the core data that represents a "node".
/// Most functional logic will not be in the struct, this is mostly a container.
#[derive(Debug)]
struct WFANode {
    /// this node's index
    node_index: usize,
    /// the sequence contained by the POA node
    sequence: Vec<u8>,
    /// contains the indices of the parent nodes, sorted
    parent_nodes: Vec<usize>
}

impl WFANode {
    /// Create a new WFANode and performs sanity checks on inputs.
    pub fn new(node_index: usize, sequence: Vec<u8>, mut parent_nodes: Vec<usize>) -> WFANode {
        parent_nodes.sort();
        WFANode {
            node_index,
            sequence,
            parent_nodes
        }
    }

    #[allow(dead_code)]
    pub fn node_index(&self) -> usize {
        self.node_index
    }

    pub fn sequence(&self) -> &[u8] {
        &self.sequence
    }

    #[allow(dead_code)]
    pub fn parent_nodes(&self) -> &[usize] {
        &self.parent_nodes
    }
}

/// Contains functionality for building a Partial-Order Alignment (POA) graph and then aligning a sequence to it.
/// Assumes that the last node added to the graph is the "target" or "destination" for mapping.
#[derive(Default)]
pub struct WFAGraph {
    /// all the nodes in the graph so far
    nodes: Vec<WFANode>,
    /// all the edges from one node to the next
    edges: Vec<Vec<usize>>
}

impl WFAGraph {
    /// Creates a new empty graph
    pub fn new() -> WFAGraph {
        WFAGraph {
            nodes: Default::default(),
            edges: Default::default()
        }
    }

    /// Constructs a WFA graph using only heterozygous variant types. This is the het-only entry point for graph construction.
    /// # Arguments
    /// * `reference` - the reference sequence that is the backbone for the graph
    /// * `variants` - a set of heterozygous variants that we are trying to assign to a read
    /// * `ref_start` - the reference start coordinate, used for offsetting variant positions; 0-based inclusive
    /// * `ref_end` - the reference end coordinate, used for offsetting variant positions; 0-based exclusive
    /// # Errors
    /// * if there are any errors from adding a node to the graph
    pub fn from_reference_variants(reference: &[u8], variants: &[Variant], ref_start: usize, ref_end: usize) -> 
        Result<(WFAGraph, NodeAlleleMap), Box<dyn std::error::Error>> {
        Self::from_reference_variants_with_hom(
            reference,
            variants,
            &[],
            ref_start,
            ref_end
        )
    }

    /// Constructs a WFA graph using both heterozygous and homozygous variant types. This is the het/hom entry point for graph construction.
    /// # Arguments
    /// * `reference` - the reference sequence that is the backbone for the graph
    /// * `variants` - a set of heterozygous variants that we are trying to assign to a read
    /// * `hom_variants` - a set of homozygous variants that we are not trying to assign to a read, but they make alignment better
    /// * `ref_start` - the reference start coordinate, used for offsetting variant positions; 0-based inclusive
    /// * `ref_end` - the reference end coordinate, used for offsetting variant positions; 0-based exclusive
    /// # Errors
    /// * if there are any errors from adding a node to the graph
    pub fn from_reference_variants_with_hom(reference: &[u8], variants: &[Variant], hom_variants: &[Variant], ref_start: usize, ref_end: usize) -> 
        Result<(WFAGraph, NodeAlleleMap), Box<dyn std::error::Error>> {
        
        let mut graph: WFAGraph = Default::default();
        let mut node_to_alleles: NodeAlleleMap = Default::default();

        let mut previous_end: usize = ref_start;
        let mut reference_index: usize;

        // this tracks nodes that need to be reconnected next, i.e. parents for the _next_ reference node
        // initial one is empty because there is no start node yet
        let mut reference_reconnect: Vec<usize> = vec![];
        // marks the alleles that should be tied to the next reference insertion
        let mut reference_alleles: Vec<(usize, u8)> = vec![];

        // this is a queue where the key is reconnect position and value is the node to reconnect at that juncture
        let mut reconnect_queue: PriorityQueue<usize, Reverse<usize>> = PriorityQueue::new();

        let mut all_variants: Vec<(&Variant, Option<usize>)> = Vec::with_capacity(variants.len() + hom_variants.len());
        for (variant_index, variant) in variants.iter().enumerate() {
            all_variants.push((variant, Some(variant_index)));
        }
        for variant in hom_variants.iter() {
            all_variants.push((variant, None));
        }
        all_variants.sort_by(|v1, v2| v1.0.position().cmp(&v2.0.position()));

        for (variant, variant_index) in all_variants.iter() {
            if variant.is_ignored() {
                // I don't think we need a trace message here for now
                continue;
            }

            // look at where this variant is
            let variant_pos: usize = variant.position() as usize;
            let ref_len: usize = variant.get_ref_len();
            if variant_pos < ref_start {
                // this variant starts before our reference block, so ignore it
                trace!("Ignoring variant starting at {} before ref_start {}", variant_pos, ref_start);
                continue;
            }
            if variant_pos + ref_len > ref_end {
                // this variant end after our reference block, so ignore it
                trace!("Ignoring variant ending at {} after ref_end {}", variant_pos+ref_len, ref_end);
                continue;
            }

            // while we have something to reconnect that reconnects BEFORE the next variant, handle it
            while (reconnect_queue.peek().unwrap_or((&usize::MAX, &Reverse(usize::MAX))).1).0 <= variant_pos {
                // get the next thing that needs to reconnect before the next variant
                let (alt_index, Reverse(alt_reconnect)) = reconnect_queue.pop().unwrap();
                assert!(alt_reconnect > previous_end);

                // first, we have to build up the reference node up until the reconnect point
                let ref_sequence: Vec<u8> = reference[previous_end..alt_reconnect].to_vec();
                reference_index = graph.add_node(ref_sequence, reference_reconnect)?;
                if !reference_alleles.is_empty() {
                    node_to_alleles.insert(reference_index, reference_alleles);
                    reference_alleles = vec![];
                }
                previous_end = alt_reconnect;
                
                // now prep the next one by marking that reference node plus the alt_index we are reconnecting
                reference_reconnect = vec![reference_index, alt_index];

                // also check if any other reconnects have an identical reconnect point
                while reconnect_queue.peek().unwrap_or((&usize::MAX, &Reverse(usize::MAX))).1.0 == alt_reconnect {
                    let (ai2, Reverse(ar2)) = reconnect_queue.pop().unwrap();
                    assert_eq!(alt_reconnect, ar2);
                    reference_reconnect.push(ai2);
                }
            }

            // at this point, any reconnections before this point have been resolved
            
            // check if the reference ended upstream of this variant (which is usually true) OR
            //   if the graph is currently empty, indicating that there is a variant at position 0 and we need a dummy start node
            if previous_end < variant_pos || graph.get_num_nodes() == 0 {
                // we need to catch up the reference, add a node representing all sequence up to this point
                let ref_sequence: Vec<u8> = reference[previous_end..variant_pos].to_vec();
                
                // now add the reference node that catches us up to this variant
                reference_index = graph.add_node(ref_sequence, reference_reconnect)?;
                if !reference_alleles.is_empty() {
                    node_to_alleles.insert(reference_index, reference_alleles);
                    reference_alleles = vec![];
                }

                // set these fields for the next reference node that gets added
                reference_reconnect = vec![reference_index];
                previous_end = variant_pos;
            } else {
                assert!(previous_end == variant_pos);
                // in this situation, we have already generated the sequence up to this variant, likely because two variants start at the same location
                // we should not have to do anything special because we already know the upstream index
            }

            // now add the alt allele(s)
            if variant.convert_index(0) != 0 {
                // allele0 is an alt, so this must be multi-allelic; basically do the same thing we would do for allele1
                // add the sequence exactly with just the immediately upstream reference node
                let alt_sequence: Vec<u8> = variant.get_truncated_allele0().to_vec();
                let parent_nodes: Vec<usize> = reference_reconnect.clone();
                let alt_index: usize = graph.add_node(alt_sequence, parent_nodes)?;
                let alt_reconnect: usize = variant_pos + ref_len;

                // also mark this alt node has having this particular allele0
                if let Some(vi) = variant_index {
                    node_to_alleles.insert(alt_index, vec![(*vi, 0)]);
                }

                // now we need to mark this new node for reconnection downstream
                reconnect_queue.push(alt_index, Reverse(alt_reconnect));
            } else {
                // the 0 allele is just reference, so add it to the reference allele set
                if let Some(vi) = variant_index {
                    reference_alleles.push((*vi, 0));
                }
            }

            // allele1 is *always* an alt, add the sequence exactly with just the immediately upstream reference node
            let alt_sequence: Vec<u8> = variant.get_truncated_allele1().to_vec();
            let parent_nodes: Vec<usize> = reference_reconnect.clone();
            let alt_index: usize = graph.add_node(alt_sequence, parent_nodes)?;
            let alt_reconnect: usize = variant_pos + ref_len;

            // also mark this alt node has having this particular allele
            if let Some(vi) = variant_index {
                node_to_alleles.insert(alt_index, vec![(*vi, 1)]);
            }

            // now we need to mark this new node for reconnection downstream
            reconnect_queue.push(alt_index, Reverse(alt_reconnect));
        }

        // reconnect everything downstream from here
        while !reconnect_queue.is_empty() {
            let (alt_index, Reverse(alt_reconnect)) = reconnect_queue.pop().unwrap();
            assert!(alt_reconnect > previous_end);
            let ref_sequence: Vec<u8> = reference[previous_end..alt_reconnect].to_vec();
            reference_index = graph.add_node(ref_sequence, reference_reconnect)?;
            if !reference_alleles.is_empty() {
                node_to_alleles.insert(reference_index, reference_alleles);
                reference_alleles = vec![];
            }
            previous_end = alt_reconnect;
            
            // now prep the next one
            reference_reconnect = vec![reference_index, alt_index];
            while reconnect_queue.peek().unwrap_or((&usize::MAX, &Reverse(usize::MAX))).1.0 == alt_reconnect {
                let (ai2, Reverse(ar2)) = reconnect_queue.pop().unwrap();
                assert_eq!(alt_reconnect, ar2);
                reference_reconnect.push(ai2);
            }
        }

        // now we just have one last reference node to add
        assert!(previous_end <= ref_end);
        let ref_sequence: Vec<u8> = reference[previous_end..ref_end].to_vec();
        graph.add_node(ref_sequence, reference_reconnect)?;

        // make sure we didn't have any loose reference alleles hanging about, I don't think this can happen unless users enter weird stuff
        assert!(reference_alleles.is_empty());

        Ok((graph, node_to_alleles))
    }

    pub fn get_num_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// Adds a node to the graph and returns its index as a Result.
    /// # Arguments
    /// * `sequence` - the vector of sequence to add with this node
    /// * `parent_nodes` - the index of any upstream nodes in the graph
    /// # Errors
    /// * if the first node inserted has parents; first node is assumed root, so this would break that assumption
    /// * if any subsequent node is parent-less; all nodes must stem from the root
    /// * if any parent node has a index >= this node's index; this is a DAG only
    pub fn add_node(&mut self, sequence: Vec<u8>, parent_nodes: Vec<usize>) -> Result<usize, Box<dyn std::error::Error>> {
        let new_index: usize = self.nodes.len();

        // sanity checks on what is being added, node wise anyways
        if new_index == 0 {
            // this is the first ndoe, it should not have any parents
            if !parent_nodes.is_empty() {
                bail!("First node must have no parent nodes.");
            }
        } else {
            // this is a non-first node, it MUST have a parent
            if parent_nodes.is_empty() {
                bail!("All nodes after the first must have at least one parent node.");
            }
            // make sure all parent node indices comes before this node
            for &pn in parent_nodes.iter() {
                if new_index <= pn {
                    bail!("All parent nodes must come before this node.");
                }
            }
        }

        // add any new edges from parents
        for &p_index in parent_nodes.iter() {
            self.edges[p_index].push(new_index);
        }
        
        // add the new node with an empty set of edges coming from it
        let new_node: WFANode = WFANode::new(new_index, sequence, parent_nodes);
        self.nodes.push(new_node);
        self.edges.push(vec![]);
        
        Ok(new_index)
    }

    /// Calculates the edit distance of `other_sequence` onto this graph, returning both the score and the traversed nodes to achieve that score.
    /// # Arguments
    /// * `other_sequence` - the sequence being aligned to this graph
    /// # Errors
    /// * if the maximum edit distance is reached; this is a safeguard from run-away loops
    pub fn edit_distance(&self, other_sequence: &[u8]) -> Result<WFAResult, Box<dyn std::error::Error>> {
        self.edit_distance_with_pruning(other_sequence, usize::MAX)
    }

    /// Calculates the edit distance of `other_sequence` onto this graph using the WFA algorithm, returning both the score and the traversed nodes to achieve that score.
    /// If multiple paths exist that are equal, all nodes along each path are returned, allowing us to mark ambiguity.
    /// This mode will prune wavefronts that fall too far behind the farthest, leading to potentially incorrect results under certain conditions.
    /// # Arguments
    /// * `other_sequence` - the sequence being aligned to this graph
    /// * `prune_distance` - if a wavefront is behind the farthest wavefront by this distance, it will be pruned; set to usize::MAX to disable pruning
    /// # Errors
    /// * if the maximum edit distance is reached; this is a safeguard from run-away loops
    pub fn edit_distance_with_pruning(&self, other_sequence: &[u8], prune_distance: usize) -> Result<WFAResult, Box<dyn std::error::Error>> {
        // We will structure the algorithm mentally such that X-axis is the graph and Y-axis is `other_sequence`.
        // This means we are iterating on columns representing characters in the graph.
        // Each column will be other_len long.

        // each node *may* have a set of active wavefronts indicating progression - HashMap of (at most) length N
        //   these active wavefronts must all be the same ED, but may not be adjacent (if one branch has a big indel) - HashMap based on start position in other_sequence
        //     each start position should have current distance into the `other_sequence` along with a set of upstream nodes - (a, Vec(b)) where a is an offset and b is index of set(s)
        let mut active_wavefronts: HashMap<usize, HashMap<isize, Vec<(usize, usize)>>> = Default::default();
        let mut next_wavefronts: HashMap<usize, HashMap<isize, Vec<(usize, usize)>>> = Default::default();
        let mut max_wavefronts: HashMap<usize, HashMap<isize, usize>> = Default::default();

        // we also need to track which nodes are traversed by our particular path
        let mut treeset_to_index: HashMap<BitVec, usize> = HashMap::default();
        let mut index_to_treeset: Vec<BitVec> = vec![];

        // we always start in node 0, so make that set
        let mut base_bitvec: BitVec = BitVec::from_elem(self.nodes.len(), false);
        base_bitvec.set(0, true);
        assert!(treeset_to_index.insert(base_bitvec.clone(), 0).is_none());
        index_to_treeset.push(base_bitvec);
        
        let base_hashset_index: usize = 0;

        // insert the starting wavefront at `other_sequence`[0] + offset = 0, with set containing just node #0
        let mut initial_wavefront: HashMap<isize, Vec<(usize, usize)>> = Default::default();
        initial_wavefront.insert(0, vec![(0, base_hashset_index)]);
        // this goes into node 0
        active_wavefronts.insert(0, initial_wavefront);

        // in a given loop, these will track any nodes that were actively moving wavefronts
        let mut min_active_wavefront: usize;
        let mut max_active_wavefront: usize;

        // each loop of WFA will increase our edit distance by 1
        let mut edit_distance: usize = 0;
        let max_edit_distance: usize = 100000;
        let mut farthest_progression: usize = 0;
        let mut min_progression: usize = 0;
        
        let mut encountered_nodes: HashMap<usize, bool> = Default::default();

        loop {
            /*
             * Outline of the core POA-WFA algorithm
             * 1. Extend all current wavefronts and create splits - this is loop 1, traversed in node order
             * 2. REMOVED  - Back-propagate maximum progressions - this is loop 2, traversed in reverse node order; 
             *      this was a net drag on the runtime AND would get wrong results on occasion
             * 3. Increase edit distance to match splits
             */
            min_active_wavefront = usize::MAX;
            max_active_wavefront = 0;

            // trace!("WFAGraph ed={} start: farthest_progression = {}, set_len = {}", edit_distance, farthest_progression, index_to_treeset.len());

            // we can iterate over our nodes in order because they are DAGs entered in order
            let mut wavefronts_scanned = 0;
            for (node_index, node) in self.nodes.iter().enumerate() {
                /* 
                 * Outline of this core extension loop:
                 * 1. Push all wavefronts for this node forward
                 * 2. Collapse all wavefronts on the same diagonal such that only the best remain.
                 *    - These get consolidated into a single WF.  If it's worse than the best so far, we remove it from consideration.
                 * 3. If a diagonal hits the end of this node's sequence, copy it into all of the children nodes.
                 *    - If it's the final node and has progress through all sequence, we instead mark it as finished.
                 * 4. Once we reach the end, generate any splits that must have ed = ed+1
                 */
                
                // check if this node has *any* active wavefronts
                if !active_wavefronts.contains_key(&node_index) {
                    continue;
                }

                // if this entry doesn't exist, set it
                if let std::collections::hash_map::Entry::Vacant(e) = encountered_nodes.entry(node_index) {
                    e.insert(true);
                    trace!("WFAGraph n#{} start: min ed = {}", node_index, edit_distance);
                }

                // update our min & max
                min_active_wavefront = min_active_wavefront.min(node_index);
                max_active_wavefront = max_active_wavefront.max(node_index);

                let node_sequence: &[u8] = node.sequence();
                let node_length: usize = node_sequence.len();
                
                // pull out the active wavefront for this node
                let mut wavefront: HashMap<isize, Vec<(usize, usize)>> = active_wavefronts.remove(&node_index).unwrap();
                let maxfront: &mut HashMap<isize, usize> = max_wavefronts.entry(node_index).or_insert_with(Default::default);

                // `other_start` represent the first position in `other_sequence` for this WF diagonal, 
                //      and it *can* be negative when we "delete" more node sequence than other sequence
                // `offset` (below) represents the offset into the current node we are comparing currently
                // if `other_start` is negative, then the corresponding `offset` values must be positive enough to overcome it (e.g. >= 0 when added)
                for (other_start, vec_waves) in wavefront.iter_mut() {
                    wavefronts_scanned += 1;

                    // first extend all wavefronts as far as possible, tracking the farthest
                    let mut max_offset: usize = 0;
                    for (offset, _hashset_index) in vec_waves.iter_mut() {
                        // get the position in `other_sequence` we are currently comparing against
                        assert!(other_start + *offset as isize >= 0);
                        let mut other_position: usize = (other_start + *offset as isize) as usize;

                        // now extend as far as we can, making sure to check for boundaries and inequality in bases
                        while *offset < node_length && 
                            other_position < other_sequence.len() && 
                            node_sequence[*offset] == other_sequence[other_position] {
                            *offset += 1;
                            other_position += 1;
                        }
                        max_offset = max_offset.max(*offset);
                    }

                    // if we go along a diagonal and max is less than we've seen before, then this is a suboptimal solution we ignore
                    let maxfront_record: &mut usize = maxfront.entry(*other_start).or_insert(0);    
                    if max_offset < *maxfront_record  || (other_start + max_offset as isize) < min_progression as isize {
                        // skips_triggered += 1;
                        // vec_waves.clear();
                        continue;
                    }
                    *maxfront_record = max_offset;

                    // double check this truth
                    assert!(other_start+max_offset as isize >= 0);
                    farthest_progression = farthest_progression.max((other_start+max_offset as isize) as usize);

                    // now collapse down everything that made it to the max_offset
                    let best_offset = max_offset;
                    let mut best_sets: Vec<usize> = vec![];
                    for &(o, s) in vec_waves.iter() {
                        if o == best_offset {
                            best_sets.push(s);
                        }
                    }

                    // remove duplicates
                    best_sets.sort();
                    best_sets.dedup();

                    let best_set: usize = if best_sets.len() > 1 {
                        // we have multiple bests, collapse them into a single index
                        let mut set_union: BitVec = BitVec::from_elem(self.nodes.len(), false);
                        for &set_index in best_sets.iter() {
                            let other_set: &BitVec = &index_to_treeset[set_index];
                            set_union.or(other_set);
                        }

                        // get the index of this entry (or create one if necessary)
                        let new_set_index: usize = match treeset_to_index.get(&set_union) {
                            Some(i) => { *i },
                            None => {
                                index_to_treeset.push(set_union.clone());
                                treeset_to_index.insert(set_union, index_to_treeset.len() - 1);
                                index_to_treeset.len() - 1
                            } 
                        };
                        new_set_index
                    } else {
                        // only one remains, just copy it
                        best_sets[0]
                    };

                    if max_offset == node_length {
                        // we are at the end of this node, do different things depending on if this is the final node or not
                        if node_index == self.nodes.len() - 1 {
                            assert!(other_start + max_offset as isize >= 0);
                            if ((other_start + max_offset as isize) as usize) < other_sequence.len() {
                                // we are *not* at the end of other sequence, but we *are* at the end of the graph
                                // now we would normally split this into three waves on this node, but only the +1 is valid in this situation
                                let node_wf: &mut HashMap<isize, Vec<(usize, usize)>> = next_wavefronts.entry(node_index).or_insert_with(Default::default);
                                // +1 on diagonal - graph does not advance, other does (other has relative insertion); other_start is one more, but offset does not increase
                                let plus_diagonal : &mut Vec<(usize, usize)> = node_wf.entry(*other_start+1).or_insert(vec![]);
                                plus_diagonal.push((max_offset, best_set));
                            } else {
                                // we are at the end of both final node and the other sequence
                                // we will handle anything below
                            }
                        } else {
                            assert!(other_start + max_offset as isize >= 0);
                            
                            // we are not in the final node, so we need to push this to successor nodes for more extension
                            // the `new_offset` tells our algorithm which base we're comparing and orients us to a diagonal
                            let new_offset: isize = other_start + max_offset as isize;
                            for &successor_index in self.edges[node_index].iter() {
                                let node_wf: &mut HashMap<isize, Vec<(usize, usize)>> = active_wavefronts.entry(successor_index).or_insert_with(Default::default);
                                let copy_diagonal: &mut Vec<(usize, usize)> = node_wf.entry(new_offset).or_insert(vec![]);
                                
                                // the successor set should include the best + the successors node index
                                let current_set: &BitVec = &index_to_treeset[best_set];
                                let mut new_set: BitVec = BitVec::from_elem(self.nodes.len(), false);
                                new_set.set(successor_index, true);
                                new_set.or(current_set);

                                // get the index of this entry (or create one if necessary)
                                let new_set_index: usize = match treeset_to_index.get(&new_set) {
                                    Some(i) => { *i },
                                    None => {
                                        index_to_treeset.push(new_set.clone());
                                        treeset_to_index.insert(new_set, index_to_treeset.len() - 1);
                                        index_to_treeset.len() - 1
                                    } 
                                };
                                copy_diagonal.push((0, new_set_index));
                            }
                        }
                    } else {
                         // now we split this into three waves on this node
                         let node_wf: &mut HashMap<isize, Vec<(usize, usize)>> = next_wavefronts.entry(node_index).or_insert_with(Default::default);
                         // -1 on diagonal - graph advances, other does not (other has relative deletion); other_start is one less, but the offset increases still
                         let minus_diagonal: &mut Vec<(usize, usize)> = node_wf.entry(*other_start-1).or_insert(vec![]);
                         minus_diagonal.push((max_offset+1, best_set));
                         
                         // these two can only happen if sequence remains in other
                         assert!(*other_start + max_offset as isize >= 0);
                         if ((*other_start + max_offset as isize) as usize) < other_sequence.len() {
                             // +0 on diagonal - both node and other advance with mismatch; other_start does not change, but offset increases +1
                             let zero_diagonal: &mut Vec<(usize, usize)> = node_wf.entry(*other_start).or_insert(vec![]);
                             zero_diagonal.push((max_offset+1, best_set));
                             
                             // +1 on diagonal - graph does not advance, other does (other has relative insertion); other_start is one more, but offset does not increase
                             let plus_diagonal : &mut Vec<(usize, usize)> = node_wf.entry(*other_start+1).or_insert(vec![]);
                             plus_diagonal.push((max_offset, best_set));
                         }
                    }
                }
                
                if node_index == self.nodes.len() - 1 {
                    // we are at the last node, check if we reached the end
                    // this will store any final results
                    let mut final_hashsets: Vec<usize> = vec![];
                    for (other_start, vec_waves) in wavefront.iter() {
                        for &(offset, hashset_index) in vec_waves.iter() {
                            // if we are at the end of the node AND our sequence
                            assert!(other_start + offset as isize >= 0);
                            if offset == node_length && (other_start + offset as isize) as usize == other_sequence.len() {
                                final_hashsets.push(hashset_index);
                            }
                        }
                    }

                    if !final_hashsets.is_empty() {
                        // we've reached the end through one or more means, collapse and return
                        // remove duplicates
                        final_hashsets.sort();
                        final_hashsets.dedup();
                    
                        // now collapse the non-duplicates if necessary
                        let best_set: usize = if final_hashsets.len() > 1 {
                            // we have more than one, merge them all together and then return that one
                            let mut set_union: BitVec = BitVec::from_elem(self.nodes.len(), false);
                            for &set_index in final_hashsets.iter() {
                                let other_set: &BitVec = &index_to_treeset[set_index];
                                set_union.or(other_set);
                            }

                            // get the index of this entry (or create one if necessary)
                            let new_set_index: usize = match treeset_to_index.get(&set_union) {
                                Some(i) => { *i },
                                None => {
                                    index_to_treeset.push(set_union.clone());
                                    treeset_to_index.insert(set_union, index_to_treeset.len() - 1);
                                    index_to_treeset.len() - 1
                                } 
                            };
                            new_set_index
                        } else {
                            // only one exists, just copy it
                            final_hashsets[0]
                        };

                        let sorted_traversed_nodes: Vec<usize> = index_to_treeset[best_set].iter()
                            .enumerate()
                            .filter(|(_i, b)| *b)
                            .map(|(i, _b)| i)
                            .collect();
                        return Ok(WFAResult {
                            score: edit_distance,
                            traversed_nodes: sorted_traversed_nodes
                        });
                    }
                }
            }

            // end of loop - increase ED and update active wavefronts
            edit_distance += 1;
            active_wavefronts = next_wavefronts;
            next_wavefronts = Default::default();

            if farthest_progression > prune_distance {
                min_progression = farthest_progression - prune_distance;
            }

            trace!("edit_distance => {}, wave_fronts scanned => {}, active_indices={}..{}", edit_distance, wavefronts_scanned, min_active_wavefront, max_active_wavefront);

            // safety while debugging
            if edit_distance > max_edit_distance {
                bail!("Max_edit_distance ({}) reached during WFA solving", max_edit_distance);
            }
        }
    }
}

/// Container for POA results.
#[derive(Debug, Eq, PartialEq)]
pub struct WFAResult {
    /// The score of the best match from the alignment
    score: usize,
    /// Nodes that were traversed to get this best match; conflicting node results indicate a tie in which branch should be traversed.
    traversed_nodes: Vec<usize>
}

impl WFAResult {
    pub fn score(&self) -> usize {
        self.score
    }
    
    pub fn traversed_nodes(&self) -> &[usize] {
        &self.traversed_nodes
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_single_node() {
        // create a new graph and add a single node to it
        let mut graph: WFAGraph = WFAGraph::new();
        let v1: Vec<u8> = vec![0, 1, 2, 4, 5];
        graph.add_node(v1.clone(), vec![]).unwrap();

        // test sequences
        let v2: Vec<u8> = vec![0, 1, 3, 4, 5];
        let v3: Vec<u8> = vec![1, 2, 3, 5];
        let v4: Vec<u8> = vec![];

        // check the nodes in the first one, but just score for the rest
        assert_eq!(graph.edit_distance(&v1).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0] });
        assert_eq!(graph.edit_distance(&v2).unwrap().score(), 1);
        assert_eq!(graph.edit_distance(&v3).unwrap().score(), 2);
        assert_eq!(graph.edit_distance(&v4).unwrap().score(), 5);
    }

    #[test]
    fn test_two_node_single_path() {
        // this is the base sequence
        let v1: Vec<u8> = vec![0, 1, 2, 4, 5];
        for split_point in 0..v1.len() {
            // split the sequence at various points and verify everything is still correct
            let mut graph: WFAGraph = WFAGraph::new();
            graph.add_node(v1[0..split_point].to_vec(), vec![]).unwrap();
            graph.add_node(v1[split_point..].to_vec(), vec![0]).unwrap();

            // test sequences
            let v2: Vec<u8> = vec![0, 1, 3, 4, 5];
            let v3: Vec<u8> = vec![1, 2, 3, 5];
            let v4: Vec<u8> = vec![];

            // check the nodes in the first one, but just score for the rest
            assert_eq!(graph.edit_distance(&v1).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 1] });
            assert_eq!(graph.edit_distance(&v2).unwrap(), WFAResult { score: 1, traversed_nodes: vec![0, 1] });
            assert_eq!(graph.edit_distance(&v3).unwrap(), WFAResult { score: 2, traversed_nodes: vec![0, 1] });
            assert_eq!(graph.edit_distance(&v4).unwrap(), WFAResult { score: 5, traversed_nodes: vec![0, 1] });
        }
    }

    #[test]
    fn test_basic_variant() {
        // create a graph where index 2 is an SNV change to either 2 or 3
        let mut graph: WFAGraph = WFAGraph::new();
        let v1: Vec<u8> = vec![0, 1, 2, 4, 5];
        graph.add_node(v1[0..2].to_vec(), vec![]).unwrap();
        graph.add_node(vec![2], vec![0]).unwrap();
        graph.add_node(vec![3], vec![0]).unwrap();
        graph.add_node(v1[3..].to_vec(), vec![1, 2]).unwrap();

        // test sequences
        // v2 is in the graph
        let v2: Vec<u8> = vec![0, 1, 3, 4, 5];
        // v3 is 2 away from v1, 3 away from v2
        let v3: Vec<u8> = vec![1, 2, 3, 5];
        // 5 away from both
        let v4: Vec<u8> = vec![];
        // 1 away from both
        let v5: Vec<u8> = vec![0, 1, 4, 5];

        // check the nodes in the first one, but just score for the rest
        assert_eq!(graph.edit_distance(&v1).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 1, 3] });
        assert_eq!(graph.edit_distance(&v2).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 2, 3] });
        assert_eq!(graph.edit_distance(&v3).unwrap(), WFAResult { score: 2, traversed_nodes: vec![0, 1, 3] });
        assert_eq!(graph.edit_distance(&v4).unwrap(), WFAResult { score: 5, traversed_nodes: vec![0, 1, 2, 3] });
        assert_eq!(graph.edit_distance(&v5).unwrap(), WFAResult { score: 1, traversed_nodes: vec![0, 1, 2, 3] });
    }

    #[test]
    fn test_triple_split() {
        let v1 = vec![0, 1, 2, 3, 4, 5];
        let v2 = vec![0, 1, 2, 4, 4, 5];
        let v3 = vec![0, 1,    4, 4, 5];

        // this construct splits the middle into 3 separate alleles (2, 3), (2, 4), and (-, 4)
        let mut graph: WFAGraph = WFAGraph::new();
        let root = graph.add_node(v1[0..2].to_vec(), vec![]).unwrap();
        let s1 = graph.add_node(v1[2..4].to_vec(), vec![root]).unwrap();
        let s2 = graph.add_node(v2[2..4].to_vec(), vec![root]).unwrap();
        let s3 = graph.add_node(v3[2..3].to_vec(), vec![root]).unwrap();
        let tail = graph.add_node(v1[4..].to_vec(), vec![s1, s2, s3]).unwrap();

        assert_eq!(graph.edit_distance(&v1).unwrap(), WFAResult { score: 0, traversed_nodes: vec![root, s1, tail] });
        assert_eq!(graph.edit_distance(&v2).unwrap(), WFAResult { score: 0, traversed_nodes: vec![root, s2, tail] });
        assert_eq!(graph.edit_distance(&v3).unwrap(), WFAResult { score: 0, traversed_nodes: vec![root, s3, tail] });
    }

    #[test]
    fn test_nested_split() {
        let v1 = vec![0, 1, 2, 3, 4, 5];
        let v2 = vec![0, 1, 2, 4, 4, 5];
        let v3 = vec![0, 1,    4, 4, 5];

        // this construct pairs2 (2, 4) with (-, 4)
        let mut graph: WFAGraph = WFAGraph::new();
        let root = graph.add_node(v1[0..2].to_vec(), vec![]).unwrap();
        // (2, 3) still alone
        let s1 = graph.add_node(v1[2..4].to_vec(), vec![root]).unwrap();
        // s2 contains just (2, )
        let s2 = graph.add_node(v2[2..3].to_vec(), vec![root]).unwrap();
        // s3 contains just (4, ), but it allows you to come directly from root to enable the deletion in v3
        let s3 = graph.add_node(v2[3..4].to_vec(), vec![root, s2]).unwrap();
        // tail is the same
        let tail = graph.add_node(v1[4..].to_vec(), vec![s1, s3]).unwrap();

        assert_eq!(graph.edit_distance(&v1).unwrap(), WFAResult { score: 0, traversed_nodes: vec![root, s1, tail] });
        assert_eq!(graph.edit_distance(&v2).unwrap(), WFAResult { score: 0, traversed_nodes: vec![root, s2, s3, tail] });
        assert_eq!(graph.edit_distance(&v3).unwrap(), WFAResult { score: 0, traversed_nodes: vec![root, s3, tail] });
    }

    #[test]
    fn test_double_split() {
        let v1 = vec![0, 1, 2, 3, 4, 5];
        let v2 = vec![0, 1, 2, 4, 4, 5];
        let v3 = vec![0, 1,    4, 4, 5];

        // this construct separate the deletion event from the SNV event with a "gap" in the middle
        let mut graph: WFAGraph = WFAGraph::new();
        let root = graph.add_node(v1[0..2].to_vec(), vec![]).unwrap();
        // s1 just contains (2, ) 
        let s1 = graph.add_node(v1[2..3].to_vec(), vec![root]).unwrap();
        // s2 is an empty join
        let s2 = graph.add_node(vec![], vec![root, s1]).unwrap();
        // s3 contains just (3, )
        let s3 = graph.add_node(v1[3..4].to_vec(), vec![s2]).unwrap();
        // s4 contains just (4, )
        let s4 = graph.add_node(v2[3..4].to_vec(), vec![s2]).unwrap();
        // tail picks up the last two bases
        let tail = graph.add_node(v1[4..].to_vec(), vec![s3, s4]).unwrap();

        assert_eq!(graph.edit_distance(&v1).unwrap(), WFAResult { score: 0, traversed_nodes: vec![root, s1, s2, s3, tail] });
        assert_eq!(graph.edit_distance(&v2).unwrap(), WFAResult { score: 0, traversed_nodes: vec![root, s1, s2, s4, tail] });
        assert_eq!(graph.edit_distance(&v3).unwrap(), WFAResult { score: 0, traversed_nodes: vec![root, s2, s4, tail] });
    }

    #[test]
    fn test_overlapping_split() {
        let v1 = vec![0, 1, 2, 3, 4, 5];
        let v2 = vec![0,       3, 4, 5];
        let v3 = vec![0, 1,       4, 5];

        /*
        Graph structure of overlapping splits, both of which delete the "2":
               ->   ->   -v
        0 -> 1 -> 2 -> 3 -> 4,5
          ->   ->   -^
        */
        let mut graph: WFAGraph = WFAGraph::new();
        let root = graph.add_node(v1[0..1].to_vec(), vec![]).unwrap();
        // represents the (1, )
        let s1 = graph.add_node(v1[1..2].to_vec(), vec![root]).unwrap();
        // represents the (2, )
        let s2 = graph.add_node(v1[2..3].to_vec(), vec![s1]).unwrap();
        // represents the (3, )
        let s3 = graph.add_node(v1[3..4].to_vec(), vec![root, s2]).unwrap();
        // represent the tail
        let tail = graph.add_node(v1[4..].to_vec(), vec![s1, s3]).unwrap();

        assert_eq!(graph.edit_distance(&v1).unwrap(), WFAResult { score: 0, traversed_nodes: vec![root, s1, s2, s3, tail] });
        assert_eq!(graph.edit_distance(&v2).unwrap(), WFAResult { score: 0, traversed_nodes: vec![root, s3, tail] });
        assert_eq!(graph.edit_distance(&v3).unwrap(), WFAResult { score: 0, traversed_nodes: vec![root, s1, tail] });
    }

    #[test]
    fn test_simple_snv() {
        let reference = "AAA".as_bytes();
        let variants = vec![Variant::new_snv(0, 1, "A".as_bytes().to_vec(), "C".as_bytes().to_vec(), 0, 1)];

        let (graph, node_to_alleles): (WFAGraph, NodeAlleleMap) = 
            WFAGraph::from_reference_variants(&reference, &variants, 0, reference.len()).unwrap();

        // check the alignments first
        assert_eq!(graph.get_num_nodes(), 4);
        // remember ALT alleles get added before reference alleles
        assert_eq!(graph.edit_distance(&reference).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 2, 3] });
        assert_eq!(graph.edit_distance("ACA".as_bytes()).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 1, 3] });
        assert_eq!(graph.edit_distance("AA".as_bytes()).unwrap(), WFAResult { score: 1, traversed_nodes: vec![0, 1, 2, 3] });

        // now check our lookup tables
        assert_eq!(*node_to_alleles.get(&0).unwrap_or(&vec![]), vec![]);
        assert_eq!(*node_to_alleles.get(&1).unwrap_or(&vec![]), vec![(0, 1)]);
        assert_eq!(*node_to_alleles.get(&2).unwrap_or(&vec![]), vec![(0, 0)]);
        assert_eq!(*node_to_alleles.get(&3).unwrap_or(&vec![]), vec![]);
    }

    #[test]
    fn test_multiple_variants() {
        // two A>C SNVs with reference in between
        let reference = "AAAAA".as_bytes();
        let variants = vec![
            // vcf_index, position, allele0, allele1, index_allele0, index_allele1
            Variant::new_snv(0, 1, "A".as_bytes().to_vec(), "C".as_bytes().to_vec(), 0, 1),
            Variant::new_snv(0, 3, "A".as_bytes().to_vec(), "C".as_bytes().to_vec(), 0, 1)
        ];

        let (graph, node_to_alleles): (WFAGraph, NodeAlleleMap) = 
            WFAGraph::from_reference_variants(&reference, &variants, 0, reference.len()).unwrap();
        
        // check the alignments first
        assert_eq!(graph.get_num_nodes(), 7);
        // remember ALT alleles get added before reference alleles
        /*
        REF: 0 -> 2 -> 3 -> 5 -> 6
        ALT:   -> 1 -^   -> 4 -^
         */
        assert_eq!(graph.edit_distance(&reference).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 2, 3, 5, 6] });
        assert_eq!(graph.edit_distance("ACAAA".as_bytes()).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 1, 3, 5, 6] });
        assert_eq!(graph.edit_distance("AAACA".as_bytes()).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 2, 3, 4, 6] });
        assert_eq!(graph.edit_distance("ACACA".as_bytes()).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 1, 3, 4, 6] });
        assert_eq!(graph.edit_distance("AAA".as_bytes()).unwrap(), WFAResult { score: 2, traversed_nodes: vec![0, 1, 2, 3, 4, 5, 6] });
        assert_eq!(graph.edit_distance("AGAGA".as_bytes()).unwrap(), WFAResult { score: 2, traversed_nodes: vec![0, 1, 2, 3, 4, 5, 6] });

        // check some mismatches on reference real quick also
        assert_eq!(graph.edit_distance("GAAAA".as_bytes()).unwrap(), WFAResult { score: 1, traversed_nodes: vec![0, 2, 3, 5, 6] });
        assert_eq!(graph.edit_distance("ACAGAA".as_bytes()).unwrap(), WFAResult { score: 1, traversed_nodes: vec![0, 1, 3, 5, 6] });

        // now check our lookup tables
        assert_eq!(*node_to_alleles.get(&0).unwrap_or(&vec![]), vec![]);
        assert_eq!(*node_to_alleles.get(&1).unwrap_or(&vec![]), vec![(0, 1)]);
        assert_eq!(*node_to_alleles.get(&2).unwrap_or(&vec![]), vec![(0, 0)]);
        assert_eq!(*node_to_alleles.get(&3).unwrap_or(&vec![]), vec![]);
        assert_eq!(*node_to_alleles.get(&4).unwrap_or(&vec![]), vec![(1, 1)]);
        assert_eq!(*node_to_alleles.get(&5).unwrap_or(&vec![]), vec![(1, 0)]);
        assert_eq!(*node_to_alleles.get(&6).unwrap_or(&vec![]), vec![]);
    }

    #[test]
    fn test_overlapping_variants() {
        let reference = "ACGTA".as_bytes();
        let variants = vec![
            // vcf_index, position, ref_len, allele0, allele1, index_allele0, index_allele1
            Variant::new_deletion(0, 1, 2, "CG".as_bytes().to_vec(), "C".as_bytes().to_vec(), 0, 1),
            Variant::new_deletion(0, 2, 2, "GT".as_bytes().to_vec(), "G".as_bytes().to_vec(), 0, 1)
        ];

        let (graph, node_to_alleles): (WFAGraph, NodeAlleleMap) = 
            WFAGraph::from_reference_variants(&reference, &variants, 0, reference.len()).unwrap();
        
        // check the alignments first
        assert_eq!(graph.get_num_nodes(), 7);
        // remember ALT alleles get added before reference alleles
        /*
        REF: 0 -> 2 -> 4 -> 5 -> 6
        ALT:   -> 1      -^
                    -> 3      -^
         */
        assert_eq!(graph.edit_distance(&reference).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 2, 4, 5, 6] });
        assert_eq!(graph.edit_distance("ACTA".as_bytes()).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 1, 5, 6] });
        assert_eq!(graph.edit_distance("ACGA".as_bytes()).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 2, 3, 6] });
        // ACGTA -> AGTA, 1 ed; ACTA -> AGTA, 1 ed; so first allele is ambiguous, other should be reference
        assert_eq!(graph.edit_distance("AGTA".as_bytes()).unwrap(), WFAResult { score: 1, traversed_nodes: vec![0, 1, 2, 4, 5, 6] });
        // AA -> AGTA, 2 ed; AA -> ACTA, 2 ed; AA -> ACGTA, 3 ed; but this should still lead to fully ambiguous solution without node 4 (it's skipped in both 2-ed paths)
        assert_eq!(graph.edit_distance("AA".as_bytes()).unwrap(), WFAResult { score: 2, traversed_nodes: vec![0, 1, 2, 3, 5, 6] });
        
        // now check our lookup tables
        assert_eq!(*node_to_alleles.get(&0).unwrap_or(&vec![]), vec![]);
        assert_eq!(*node_to_alleles.get(&1).unwrap_or(&vec![]), vec![(0, 1)]);
        assert_eq!(*node_to_alleles.get(&2).unwrap_or(&vec![]), vec![(0, 0)]);
        assert_eq!(*node_to_alleles.get(&3).unwrap_or(&vec![]), vec![(1, 1)]);
        assert_eq!(*node_to_alleles.get(&4).unwrap_or(&vec![]), vec![(1, 0)]);
        assert_eq!(*node_to_alleles.get(&5).unwrap_or(&vec![]), vec![]);
        assert_eq!(*node_to_alleles.get(&6).unwrap_or(&vec![]), vec![]);
    }

    #[test]
    fn test_identical_insertions() {
        // this can happen if we pair something like DV with pbmm2 and there is a "large" insertion
        let reference = "ACGTA".as_bytes();
        let variants = vec![
            // vcf_index, position, allele0, allele1, index_allele0, index_allele1
            Variant::new_insertion(0, 2, "G".as_bytes().to_vec(), "GT".as_bytes().to_vec(), 0, 1),
            Variant::new_insertion(1, 2, "G".as_bytes().to_vec(), "GT".as_bytes().to_vec(), 0, 1)
        ];

        let (graph, node_to_alleles): (WFAGraph, NodeAlleleMap) = 
            WFAGraph::from_reference_variants(&reference, &variants, 0, reference.len()).unwrap();
        
        // check the alignments first
        assert_eq!(graph.get_num_nodes(), 5);
        /* 
        REF: 0 -> 3 -> 4
        ALT:   -> 1 -^
               -> 2 -^
        */
        assert_eq!(graph.edit_distance(&reference).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 3, 4] });
        assert_eq!(graph.edit_distance("ACGTTA".as_bytes()).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 1, 2, 4] });
        // tests an insertion of the wrong character, this should be fully ambiguous
        assert_eq!(graph.edit_distance("ACGATA".as_bytes()).unwrap(), WFAResult { score: 1, traversed_nodes: vec![0, 1, 2, 3, 4] });
        
        // now check our lookup tables
        assert_eq!(*node_to_alleles.get(&0).unwrap_or(&vec![]), vec![]);
        assert_eq!(*node_to_alleles.get(&1).unwrap_or(&vec![]), vec![(0, 1)]);
        assert_eq!(*node_to_alleles.get(&2).unwrap_or(&vec![]), vec![(1, 1)]);
        assert_eq!(*node_to_alleles.get(&3).unwrap_or(&vec![]), vec![(0, 0), (1, 0)]); // this one has BOTH reference alleles
        assert_eq!(*node_to_alleles.get(&4).unwrap_or(&vec![]), vec![]);
    }

    #[test]
    fn test_multiallelic_indel() {
        // tests GT -> G/GTT
        let reference = "ACGTA".as_bytes();
        let variants = vec![
            // vcf_index, position, ref_len, allele0, allele1, index_allele0, index_allele1
            Variant::new_indel(0, 2, 2, "G".as_bytes().to_vec(), "GTT".as_bytes().to_vec(), 1, 2)
        ];

        let (graph, node_to_alleles): (WFAGraph, NodeAlleleMap) = 
            WFAGraph::from_reference_variants(&reference, &variants, 0, reference.len()).unwrap();
        
        // check the alignments first
        assert_eq!(graph.get_num_nodes(), 5);
        /* 
        REF: 0 -> 3 -> 4
        ALT:   -> 1 -^
               -> 2 -^
        */
        assert_eq!(graph.edit_distance(&reference).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 3, 4] });
        assert_eq!(graph.edit_distance("ACGA".as_bytes()).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 1, 4] });
        assert_eq!(graph.edit_distance("ACGTTA".as_bytes()).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 2, 4] });
        // we can't really get an ambiguous allele assignment here, but we can do ambiguous branches
        assert_eq!(graph.edit_distance("ACGGA".as_bytes()).unwrap(), WFAResult { score: 1, traversed_nodes: vec![0, 1, 3, 4] });
        assert_eq!(graph.edit_distance("ACGGTA".as_bytes()).unwrap(), WFAResult { score: 1, traversed_nodes: vec![0, 2, 3, 4] });
        
        // now check our lookup tables
        assert_eq!(*node_to_alleles.get(&0).unwrap_or(&vec![]), vec![]);
        assert_eq!(*node_to_alleles.get(&1).unwrap_or(&vec![]), vec![(0, 0)]); // ALT0
        assert_eq!(*node_to_alleles.get(&2).unwrap_or(&vec![]), vec![(0, 1)]); // ALT1
        assert_eq!(*node_to_alleles.get(&3).unwrap_or(&vec![]), vec![]); // we still need to inject a reference allele, but it shouldn't have an assoc. with the variant
        assert_eq!(*node_to_alleles.get(&4).unwrap_or(&vec![]), vec![]);
    }

    #[test]
    fn test_partial_reference() {
        // we prepended and appended "AA" to the basic SNV test
        let reference = "AAAAAAA".as_bytes();
        // variant coordinate shifted +2
        let variants = vec![Variant::new_snv(0, 3, "A".as_bytes().to_vec(), "C".as_bytes().to_vec(), 0, 1)];

        let (graph, node_to_alleles): (WFAGraph, NodeAlleleMap) = 
            WFAGraph::from_reference_variants(&reference, &variants, 2, reference.len()-2).unwrap();

        // check the alignments first
        assert_eq!(graph.get_num_nodes(), 4);
        // remember ALT alleles get added before reference alleles
        assert_eq!(graph.edit_distance(&reference[2..5]).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 2, 3] });
        assert_eq!(graph.edit_distance("ACA".as_bytes()).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 1, 3] });
        assert_eq!(graph.edit_distance("AA".as_bytes()).unwrap(), WFAResult { score: 1, traversed_nodes: vec![0, 1, 2, 3] });

        // now check our lookup tables
        assert_eq!(*node_to_alleles.get(&0).unwrap_or(&vec![]), vec![]);
        assert_eq!(*node_to_alleles.get(&1).unwrap_or(&vec![]), vec![(0, 1)]);
        assert_eq!(*node_to_alleles.get(&2).unwrap_or(&vec![]), vec![(0, 0)]);
        assert_eq!(*node_to_alleles.get(&3).unwrap_or(&vec![]), vec![]);
    }

    #[test]
    fn test_complex_problem() {
        let reference = "AACGTTGACGTCC".as_bytes(); // we skip 2 on the front, 1 on the back
        let variants = vec![
            // vcf_index, position, ref_len, allele0, allele1, index_allele0, index_allele1
            // 3: GTTG>G
            Variant::new_deletion(0, 3, 4, "GTTG".as_bytes().to_vec(), "G".as_bytes().to_vec(), 0, 1),
            // 4: TT>T
            Variant::new_deletion(0, 4, 2, "TT".as_bytes().to_vec(), "T".as_bytes().to_vec(), 0, 1),
            // vcf_index, position, allele0, allele1, index_allele0, index_allele1
            // 6: G>[A,C]
            Variant::new_snv(0, 6, "A".as_bytes().to_vec(), "C".as_bytes().to_vec(), 1, 2)
        ];

        let (graph, node_to_alleles): (WFAGraph, NodeAlleleMap) = 
            WFAGraph::from_reference_variants(&reference, &variants, 2, reference.len()-1).unwrap();

        // check the alignments first
        assert_eq!(graph.get_num_nodes(), 9);
        /* 
        // this tests the scenario where one variant rejoins at the start of another
        // here, the TT>T joins right before the G>[A,C] so node 3 and 4 are *both* parents of 5, 6, and 7
        REF: 0 -> 2 -> 4 -> 7 -> 8
        ALT:   -> 1           -^
                    -> 3 -|
                         -> 5 -^
                         -> 6 -^
        */
        // remember ALT alleles get added before reference alleles
        assert_eq!(graph.edit_distance("CGTTGACGTC".as_bytes()).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 2, 4, 7, 8] }); //reference
        assert_eq!(graph.edit_distance("CGACGTC".as_bytes()).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 1, 8] }); //first variant only
        assert_eq!(graph.edit_distance("CGTGACGTC".as_bytes()).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 2, 3, 7, 8] }); //second variant only
        assert_eq!(graph.edit_distance("CGTTAACGTC".as_bytes()).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 2, 4, 5, 8] }); //third-0 variant only
        assert_eq!(graph.edit_distance("CGTTCACGTC".as_bytes()).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 2, 4, 6, 8] }); //third-1 variant only
        assert_eq!(graph.edit_distance("CGTAACGTC".as_bytes()).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 2, 3, 5, 8] }); //second and third-0 variants
        assert_eq!(graph.edit_distance("CGTCACGTC".as_bytes()).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 2, 3, 6, 8] }); //second and third-1 variants

        // now try some inexacts
        assert_eq!(graph.edit_distance("CGGACGTC".as_bytes()).unwrap(), WFAResult { score: 1, traversed_nodes: vec![0, 1, 2, 3, 7, 8] }); //delete both Ts, ambiguous deletion call
        assert_eq!(graph.edit_distance("CGTACGTC".as_bytes()).unwrap(), WFAResult { score: 1, traversed_nodes: vec![0, 1, 2, 3, 5, 6, 7, 8] }); //delete TG, ambiguous in a lot of ways

        // now check our lookup tables
        assert_eq!(*node_to_alleles.get(&0).unwrap_or(&vec![]), vec![]);
        assert_eq!(*node_to_alleles.get(&1).unwrap_or(&vec![]), vec![(0, 1)]);
        assert_eq!(*node_to_alleles.get(&2).unwrap_or(&vec![]), vec![(0, 0)]);
        assert_eq!(*node_to_alleles.get(&3).unwrap_or(&vec![]), vec![(1, 1)]);
        assert_eq!(*node_to_alleles.get(&4).unwrap_or(&vec![]), vec![(1, 0)]);
        assert_eq!(*node_to_alleles.get(&5).unwrap_or(&vec![]), vec![(2, 0)]);
        assert_eq!(*node_to_alleles.get(&6).unwrap_or(&vec![]), vec![(2, 1)]);
        assert_eq!(*node_to_alleles.get(&7).unwrap_or(&vec![]), vec![]);
        assert_eq!(*node_to_alleles.get(&8).unwrap_or(&vec![]), vec![]);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // After here are mostly edge case bug tests
    ////////////////////////////////////////////////////////////////////////////////
    
    // tests when a variant starts before the provided reference, we should ignore it basically
    #[test]
    fn test_variant_before_start() {
        // the first 10 bases are ignored here "NNNNNNNNNA"
        let reference = "NNNNNNNNNAACGTA".as_bytes();
        let ref_start: usize = 10;
        // the first one should be ignored, the second one included
        let variants = vec![
            // vcf_index, position, allele0, allele1, index_allele0, index_allele1
            Variant::new_snv(0, (ref_start-1) as i64, "A".as_bytes().to_vec(), "T".as_bytes().to_vec(), 0, 1),
            Variant::new_snv(0, ref_start as i64, "A".as_bytes().to_vec(), "T".as_bytes().to_vec(), 0, 1)
        ];

        let (graph, node_to_alleles): (WFAGraph, NodeAlleleMap) = 
            WFAGraph::from_reference_variants(&reference, &variants, ref_start, reference.len()).unwrap();
        
        // check the alignments first
        assert_eq!(graph.get_num_nodes(), 4);

        // now check our lookup tables
        assert_eq!(*node_to_alleles.get(&0).unwrap_or(&vec![]), vec![]); // first node is empty
        assert_eq!(*node_to_alleles.get(&1).unwrap_or(&vec![]), vec![(1, 1)]); // second is the alternate for the SECOND variant (first variant is ignored)
        assert_eq!(*node_to_alleles.get(&2).unwrap_or(&vec![]), vec![(1, 0)]); // third is the reference for the SECOND variant
        assert_eq!(*node_to_alleles.get(&3).unwrap_or(&vec![]), vec![]); // last is just more reference
    }

    // tests when a variant goes past the provided reference, we should ignore it basically
    #[test]
    fn test_span_ref_end() {
        // tests GT -> G/GTT
        let reference = "ACGTA".as_bytes();
        let variants = vec![
            // vcf_index, position, ref_len, allele0, allele1, index_allele0, index_allele1
            Variant::new_deletion(0, 3, 3, "TAG".as_bytes().to_vec(), "T".as_bytes().to_vec(), 0, 1)
        ];

        let (graph, node_to_alleles): (WFAGraph, NodeAlleleMap) = 
            WFAGraph::from_reference_variants(&reference, &variants, 0, reference.len()).unwrap();
        
        // check the alignments first
        assert_eq!(graph.get_num_nodes(), 1);

        // now check our lookup tables
        assert_eq!(*node_to_alleles.get(&0).unwrap_or(&vec![]), vec![]);
    }

    #[test]
    fn test_hom_variants() {
        // add a hom variant at base 1; these can be traversed, but provide no index lookup because they are not a "variant"
        let reference = "AAAAA".as_bytes();
        let variants = vec![Variant::new_snv(0, 3, "A".as_bytes().to_vec(), "C".as_bytes().to_vec(), 0, 1)];
        let hom_variants = vec![Variant::new_snv(0, 1, "A".as_bytes().to_vec(), "C".as_bytes().to_vec(), 0, 1)];

        let (graph, node_to_alleles): (WFAGraph, NodeAlleleMap) = 
            WFAGraph::from_reference_variants_with_hom(&reference, &variants, &hom_variants, 0, reference.len()).unwrap();

        // check the alignments first
        assert_eq!(graph.get_num_nodes(), 7);
        // remember ALT alleles get added before reference alleles
        assert_eq!(graph.edit_distance("AAAAA".as_bytes()).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 2, 3, 5, 6] });
        assert_eq!(graph.edit_distance("ACAAA".as_bytes()).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 1, 3, 5, 6] });
        assert_eq!(graph.edit_distance("ACACA".as_bytes()).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 1, 3, 4, 6] });
        assert_eq!(graph.edit_distance("ACAA".as_bytes()).unwrap(), WFAResult { score: 1, traversed_nodes: vec![0, 1, 3, 4, 5, 6] });

        // now check our lookup tables
        assert_eq!(*node_to_alleles.get(&0).unwrap_or(&vec![]), vec![]);
        assert_eq!(*node_to_alleles.get(&1).unwrap_or(&vec![]), vec![]);
        assert_eq!(*node_to_alleles.get(&2).unwrap_or(&vec![]), vec![]);
        assert_eq!(*node_to_alleles.get(&3).unwrap_or(&vec![]), vec![]);
        assert_eq!(*node_to_alleles.get(&4).unwrap_or(&vec![]), vec![(0, 1)]);
        assert_eq!(*node_to_alleles.get(&5).unwrap_or(&vec![]), vec![(0, 0)]);
        assert_eq!(*node_to_alleles.get(&6).unwrap_or(&vec![]), vec![]);
    }

    #[test]
    fn test_variant_at_start() {
        let reference = "AAA".as_bytes();
        let variants = vec![Variant::new_snv(0, 0, "A".as_bytes().to_vec(), "C".as_bytes().to_vec(), 0, 1)];

        let (graph, node_to_alleles): (WFAGraph, NodeAlleleMap) = 
            WFAGraph::from_reference_variants(&reference, &variants, 0, reference.len()).unwrap();

        // check the alignments first
        assert_eq!(graph.get_num_nodes(), 4);
        // remember ALT alleles get added before reference alleles
        assert_eq!(graph.edit_distance(&reference).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 2, 3] });
        assert_eq!(graph.edit_distance("CAA".as_bytes()).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 1, 3] });
        assert_eq!(graph.edit_distance("AA".as_bytes()).unwrap(), WFAResult { score: 1, traversed_nodes: vec![0, 1, 2, 3] });

        // now check our lookup tables
        assert_eq!(*node_to_alleles.get(&0).unwrap_or(&vec![]), vec![]);
        assert_eq!(*node_to_alleles.get(&1).unwrap_or(&vec![]), vec![(0, 1)]);
        assert_eq!(*node_to_alleles.get(&2).unwrap_or(&vec![]), vec![(0, 0)]);
        assert_eq!(*node_to_alleles.get(&3).unwrap_or(&vec![]), vec![]);
    }

    #[test]
    fn test_variant_at_end() {
        let reference = "AAA".as_bytes();
        let variants = vec![Variant::new_snv(0, 2, "A".as_bytes().to_vec(), "C".as_bytes().to_vec(), 0, 1)];

        let (graph, node_to_alleles): (WFAGraph, NodeAlleleMap) = 
            WFAGraph::from_reference_variants(&reference, &variants, 0, reference.len()).unwrap();

        // check the alignments first
        assert_eq!(graph.get_num_nodes(), 4);
        // remember ALT alleles get added before reference alleles
        assert_eq!(graph.edit_distance(&reference).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 2, 3] });
        assert_eq!(graph.edit_distance("AAC".as_bytes()).unwrap(), WFAResult { score: 0, traversed_nodes: vec![0, 1, 3] });
        assert_eq!(graph.edit_distance("AA".as_bytes()).unwrap(), WFAResult { score: 1, traversed_nodes: vec![0, 1, 2, 3] });

        // now check our lookup tables
        assert_eq!(*node_to_alleles.get(&0).unwrap_or(&vec![]), vec![]);
        assert_eq!(*node_to_alleles.get(&1).unwrap_or(&vec![]), vec![(0, 1)]);
        assert_eq!(*node_to_alleles.get(&2).unwrap_or(&vec![]), vec![(0, 0)]);
        assert_eq!(*node_to_alleles.get(&3).unwrap_or(&vec![]), vec![]);
    }
}