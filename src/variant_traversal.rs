
use bio::data_structures::interval_tree::IntervalTree;
use log::{debug, trace};
use simple_error::{SimpleError, bail};
use rustc_hash::FxHashSet as HashSet;

use crate::astar_phaser::AstarResult;
use crate::data_types::read_segments::{AlleleType, ReadSegment};
use crate::data_types::variants::Variant;

/// This struct contains the information necessary to solve a re-ordered variant problem.
/// It can also revert the solution back to the original variant order, which we must do at the end.
pub struct OptimizedProblem {
    /// The lookup for reordering
    index_order: Vec<usize>,
    /// The reverse lookup for reordering
    reverse_lookup: Vec<usize>,
    /// The optimized variant order
    variants: Vec<Variant>,
    /// The optimized read segments
    read_segments: IntervalTree<usize, ReadSegment>,
}

impl OptimizedProblem {
    /// Constructs a new optimized problem
    /// # Arguments
    /// * `index_order` - the lookup for reordering
    /// * `variants` - the variants in the optimized order
    /// * `read_segments` - the read segments in the optimized order
    /// # Errors
    /// * if the index order is not valid, this can be caused by duplicate indices or missing indices
    pub fn new(index_order: Vec<usize>, variants: Vec<Variant>, read_segments: IntervalTree<usize, ReadSegment>) -> Result<OptimizedProblem, SimpleError> {
        // create the reverse lookup
        let mut reverse_lookup: Vec<Option<usize>> = vec![None; variants.len()];
        for (i, index) in index_order.iter().enumerate() {
            if reverse_lookup[*index].is_some() {
                bail!("Index {index} is already in the reverse lookup");
            }
            reverse_lookup[*index] = Some(i);
        }

        // validate the reverse lookup
        for (i, &index) in reverse_lookup.iter().enumerate() {
            if index.is_none() {
                bail!("Index {} is not in the reverse lookup", i);
            }
        }

        // unwrap it into a vector
        let reverse_lookup: Vec<usize> = reverse_lookup.into_iter().map(|i| i.unwrap()).collect();

        Ok(OptimizedProblem {
            index_order,
            reverse_lookup,
            variants,
            read_segments
        })
    }

    /// Translates an A* solution to match the original variant order
    /// # Arguments
    /// * `astar_result` - the A* solution to translate, which is in the OptimizedProblem order
    /// # Returns
    /// The translated A* solution, which is in the original variant order
    pub fn translate_solution(&self, astar_result: &AstarResult) -> AstarResult {
        // reorder the astar result to match the original variant order
        let reordered_haplotype_1: Vec<AlleleType> = self.reverse_lookup.iter()
            .map(|i| astar_result.haplotype_1[*i]).collect();
        let reordered_haplotype_2: Vec<AlleleType> = self.reverse_lookup.iter()
            .map(|i| astar_result.haplotype_2[*i]).collect();

        // copy the statistics
        let reordered_statistics = astar_result.statistics.clone();
        AstarResult {
            haplotype_1: reordered_haplotype_1,
            haplotype_2: reordered_haplotype_2,
            statistics: reordered_statistics
        }
    }

    // getters
    pub fn index_order(&self) -> &[usize] {
        &self.index_order
    }

    pub fn variants(&self) -> &[Variant] {
        &self.variants
    }

    pub fn read_segments(&self) -> &IntervalTree<usize, ReadSegment> {
        &self.read_segments
    }
}

/// Optimizes the variant order for the given read segments
/// # Arguments
/// * `variants` - the variants to optimize the order for
/// * `read_segments` - the read segments to optimize the order for
/// # Returns
/// The optimized problem
/// # Errors
/// * if the variant order is not valid, this can be caused by duplicate indices or missing indices
pub fn optimize_variant_order(
    variants: &[Variant], read_segments: &IntervalTree<usize, ReadSegment>
) -> Result<OptimizedProblem, SimpleError> {
    debug!("Optimizing variant order...");
    // first, compute the overlap counts for each variant pair
    let num_variants = variants.len();
    let overlap_counts = compute_overlap_counts(num_variants, read_segments);

    // find the variant with the most overlap; this will be our starting point
    let mut max_overlap = 0;
    let mut max_index = 0;
    for (i, row) in overlap_counts.iter().enumerate() {
        // the i-th variant has a read count of row[i]
        let read_count = row[i];
        if read_count > max_overlap {
            max_overlap = read_count;
            max_index = i;
        }
    }

    // TODO: there is a marginal benefit to starting at the highest coverage variant in RNA (above), but...
    //       starting at the first variant encourages a more traditional linear ordering of variants,
    //       which may encourage a more linear traversal, co-localizing nearby variants
    //       maybe we can revisit in the future if we have something that can better enforce a pseudo-linear ordering while maintaining the performance

    // always start with the first variant
    // let max_index = 0;
    // let max_overlap = overlap_counts[0][0];

    // now we have the variant with the most overlap; we will use this as our starting point
    trace!("\tFirst index {max_index} => {max_overlap} reads");
    let mut index_order = vec![max_index];
    let mut added_variants: HashSet<usize> = index_order.iter().cloned().collect();
    let mut current_overlaps = vec![0; num_variants];

    while index_order.len() < num_variants {
        // add the counts from the most recent added variant to the current overlaps
        let last_added = *index_order.last().unwrap();
        for i in 0..num_variants {
            current_overlaps[i] += overlap_counts[last_added][i];
        }

        // find the variant with the most overlap
        let mut max_overlap = 0;
        let mut max_index = None;
        for (i, &count) in current_overlaps.iter().enumerate() {
            // skip anything that has already been added
            if added_variants.contains(&i) {
                continue;
            }

            if max_index.is_none() || count > max_overlap {
                max_overlap = count;
                max_index = Some(i);
            }
        }

        // we should always get something
        let max_index = max_index.expect("No max index found, this should not happen");
        trace!("\tAdding index {max_index} => {max_overlap} pairs");
        index_order.push(max_index);
        added_variants.insert(max_index);
    }

    // fix it to non-mutable
    let index_order = index_order;
    debug!("Optimized variant order: {:?}", index_order);

    // now we need to reorder the variants
    let reordered_variants: Vec<Variant> = index_order.iter().map(|i| variants[*i].clone()).collect();

    // now we need to reorder the read segments
    let mut reordered_read_segments: IntervalTree<usize, ReadSegment> = IntervalTree::new();
    for rs_entry in read_segments.find(0..usize::MAX) {
        // get the read segment
        let rs = rs_entry.data();

        // create the new read segment components, which needs to be re-ordered to match the variant order
        let new_alleles: Vec<AlleleType> = index_order.iter().map(|i| rs.allele(*i)).collect();
        let new_quals: Vec<u8> = index_order.iter().map(|i| rs.qual(*i)).collect();

        // build the new read segment
        let new_read_segment = ReadSegment::new(
            rs.read_name().to_string(),
            new_alleles,
            new_quals
        );

        // add the new read segment to the interval tree
        let new_range = rs.region();
        reordered_read_segments.insert(new_range, new_read_segment);
    }

    // finally, but it all together into a new optimized problem
    OptimizedProblem::new(
        index_order,
        reordered_variants,
        reordered_read_segments
    )
}

/// Computes the overlap counts for each variant pair in the read segments
/// # Arguments
/// * `read_segments` - the read segments to compute the overlap counts for
/// # Returns
/// A vector of vectors, where the i-th vector contains the overlap counts for the i-th variant.
/// Index [i][i] is the total overlap count for the i-th variant.
/// Index [i][j] is the overlap count for the i-th variant and the j-th variant.
fn compute_overlap_counts(num_variants: usize, read_segments: &IntervalTree<usize, ReadSegment>) -> Vec<Vec<u64>> {
    // initialize the overlap counts to 0; this is a 2D matrix
    let mut overlap_counts: Vec<Vec<u64>> = vec![vec![0; num_variants]; num_variants];
    for rs_entry in read_segments.find(0..usize::MAX) {
        // get the read segment and its range
        let rs = rs_entry.data();
        let rs_range = rs.region();

        // iterate over the alleles in the read segment, we only care about the set alleles
        #[allow(clippy::needless_range_loop)] // easier to read this way
        for i in rs_range.start..rs_range.end {
            if rs.allele(i).is_set() {
                // the "i" allele is set, so increment the [i][i] pair as a total count
                overlap_counts[i][i] += 1;

                // iterate over the other alleles in the read segment, only check those after this point
                #[allow(clippy::needless_range_loop)] // easier to read this way
                for j in (i+1)..rs_range.end {
                    if rs.allele(j).is_set() {
                        // add forward and backward counts for faster iteration
                        overlap_counts[i][j] += 1;
                        overlap_counts[j][i] += 1;
                    }
                }
            }
        }
    }

    // return the overlap counts
    overlap_counts
}

#[cfg(test)]
mod tests {
    use crate::writers::phase_stats::PhaseStats;

    use super::*;

    #[test]
    fn test_compute_overlap_counts() {
        // build some reads
        let mut read_segments = IntervalTree::new();
        read_segments.insert(0..3, ReadSegment::new("r1".to_string(), vec![AlleleType::Reference, AlleleType::Reference, AlleleType::Reference], vec![1, 1, 1]));
        read_segments.insert(0..3, ReadSegment::new("r2".to_string(), vec![AlleleType::Reference, AlleleType::Ambiguous, AlleleType::Reference], vec![1, 0, 1]));
        read_segments.insert(0..2, ReadSegment::new("r3".to_string(), vec![AlleleType::Alternate, AlleleType::Alternate, AlleleType::NoOverlap], vec![1, 1, 0]));

        // compute the overlap counts and check them
        let overlap_counts = compute_overlap_counts(3, &read_segments);
        assert_eq!(overlap_counts, vec![
            vec![3, 2, 2],
            vec![2, 2, 1],
            vec![2, 1, 2]
        ]);
    }

    #[test]
    fn test_optimize_variant_order() {
        // build some variants
        let variants = vec![
            Variant::new_snv(0, 100, b"A".to_vec(), b"T".to_vec(), 0, 1).unwrap(),
            Variant::new_snv(0, 200, b"G".to_vec(), b"C".to_vec(), 0, 1).unwrap(),
            Variant::new_snv(0, 300, b"T".to_vec(), b"A".to_vec(), 0, 1).unwrap(),
        ];

        // build some reads
        let mut read_segments = IntervalTree::new();
        read_segments.insert(0..3, ReadSegment::new("r1".to_string(), vec![AlleleType::Reference, AlleleType::Ambiguous, AlleleType::Reference], vec![1, 0, 1]));
        read_segments.insert(0..3, ReadSegment::new("r2".to_string(), vec![AlleleType::Reference, AlleleType::Ambiguous, AlleleType::Reference], vec![1, 0, 1]));
        read_segments.insert(1..3, ReadSegment::new("r3".to_string(), vec![AlleleType::NoOverlap, AlleleType::Alternate, AlleleType::Reference], vec![0, 1, 1]));

        // optimize the variant order and check it
        let optimized_problem = optimize_variant_order(&variants, &read_segments).unwrap();
        let expected_index_order = vec![2, 0, 1];
        let reverse_lookup = vec![1, 2, 0];
        assert_eq!(optimized_problem.index_order(), &expected_index_order);

        // check the reordered variants
        let expected_variants: Vec<Variant> = expected_index_order.iter().map(|i| variants[*i].clone()).collect();
        assert_eq!(optimized_problem.variants(), &expected_variants);

        // create our expected reordered read segments and compare
        let mut reordered_read_segments: IntervalTree<usize, ReadSegment> = IntervalTree::new();
        for rs_entry in read_segments.find(0..usize::MAX) {
            let rs = rs_entry.data();
            let new_alleles: Vec<AlleleType> = expected_index_order.iter().map(|i| rs.allele(*i)).collect();
            let new_quals: Vec<u8> = expected_index_order.iter().map(|i| rs.qual(*i)).collect();
            let new_read_segment = ReadSegment::new(rs.read_name().to_string(), new_alleles, new_quals);
            reordered_read_segments.insert(rs.region(), new_read_segment);
        }
        assert_eq!(optimized_problem.read_segments(), &reordered_read_segments);

        let astar_result = AstarResult {
            haplotype_1: vec![AlleleType::Reference, AlleleType::Reference, AlleleType::Alternate],
            haplotype_2: vec![AlleleType::Alternate, AlleleType::Alternate, AlleleType::Reference],
            statistics: PhaseStats::astar_new(0, 0, 0, 3, 3, 0, 0),
        };

        // we need to translate these to get the expected haplotypes
        let expected_h1: Vec<AlleleType> = reverse_lookup.iter().map(|&i| astar_result.haplotype_1[i].clone()).collect();
        let expected_h2: Vec<AlleleType> = reverse_lookup.iter().map(|&i| astar_result.haplotype_2[i].clone()).collect();

        let translated_astar_result = optimized_problem.translate_solution(&astar_result);
        assert_eq!(translated_astar_result.haplotype_1, expected_h1);
        assert_eq!(translated_astar_result.haplotype_2, expected_h2);
        assert_eq!(translated_astar_result.statistics, astar_result.statistics);
    }
}