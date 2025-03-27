
use std::ops::Range;

#[repr(u8)]
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq, PartialOrd, strum_macros::FromRepr)]
pub enum AlleleType {
    Reference=0,
    Alternate=1,
    Ambiguous=2,
    NoOverlap=3
}

/*
impl From<u8> for AlleleType {
    fn from(value: u8) -> Self {
        match value {
            0 => AlleleType::Reference,
            1 => AlleleType::Alternate,
            2 => AlleleType::Ambiguous,
            3 => AlleleType::NoOverlap,
            _ => AlleleType::NoOverlap
        }
    }
}
*/

/// Container for a read segment that has been converted into a variant representation
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReadSegment {
    /// the read name
    read_name: String,
    /// the actual alleles for the read, should always be 0, 1, 2 (ambiguous), or 3 (non-overlapping/undefined).
    /// anything other than 0 or 1 are basically ignored in all functions
    alleles: Vec<AlleleType>,
    /// the associated quality values for converting 0 <--> 1; undefined alleles should have qual == 0
    quals: Vec<u8>,
    /// the range that is stored in alleles and quals
    region: Range<usize>
}

impl ReadSegment {
    /// Creates a new read segment from a set of alleles and quality values.
    /// Note that the segment should have the same length as the phase block, even if the segment does not actually span the full block.
    /// # Arguments
    /// * `alleles` - the alleles, should all be 0, 1, 2 (ambiguous) or 3 (non-overlapping/undefined); the first and last defined values determine the length of the segment
    /// * `quals` - cost to convert an allele from 0 <--> 1, any undefined/ambiguous alleles should be qual = 0
    /// # Panics
    /// * if `allele.len() != quals.len()`
    pub fn new(read_name: String, alleles: Vec<AlleleType>, quals: Vec<u8>) -> ReadSegment {
        assert_eq!(alleles.len(), quals.len());

        // this used to be defined by < Ambiguous, but I think we want to keep ambiguous around now
        let (first_allele, _) = alleles.iter().enumerate()
            .find(|(_i, &a)| a < AlleleType::Ambiguous).unwrap_or((alleles.len(), &AlleleType::NoOverlap));
        let last_allele = alleles.iter().enumerate().rev()
            .find(|(_i, &a)| a < AlleleType::Ambiguous) // find the last index that is set
            .map(|(i, _a)| i+1) // increment it by one since we are using exclusive-ended Range
            .unwrap_or(alleles.len()); // if we do not find one, set it to the length

        // compile the region and cut it out of the provided data
        let region = first_allele..last_allele;
        let alleles = alleles[region.clone()].to_vec();
        let quals = quals[region.clone()].to_vec();

        ReadSegment {
            read_name,
            alleles,
            quals,
            region
        }
    }

    /// Given a collection of read segments, this will collapse them into a single one.
    /// Any ambiguous/undefined alleles will have their quality set to 0.
    /// # Arguments
    /// * `read_segments` - the reads to collapse together
    /// # Panics
    /// * if `read_segments` is empty
    /// * if `read_segments` are not all of equal length
    pub fn collapse(read_segments: &[ReadSegment]) -> ReadSegment {
        // short circuit
        assert!(!read_segments.is_empty());
        if read_segments.len() == 1 {
            return read_segments[0].clone();
        }

        // figure out how large this region can be
        let min_start = read_segments.iter().map(|rs| rs.region().start).min().unwrap();
        let max_end = read_segments.iter().map(|rs| rs.region().end).max().unwrap();

        // max_end is exclusive
        let read_name: String = read_segments[0].read_name().to_string();
        let mut alleles: Vec<AlleleType> = vec![AlleleType::NoOverlap; max_end];
        let mut quals: Vec<u8> = vec![0; max_end];

        for rs in read_segments.iter() {
            assert_eq!(read_name, rs.read_name());

            for i in min_start..max_end {
                // now get the allele and quality
                let rsa = rs.allele(i);
                let rsq = rs.qual(i);

                // if rsa is unset, we skip everything
                if rsa != AlleleType::NoOverlap {
                    if alleles[i] == AlleleType::NoOverlap {
                        alleles[i] = rsa;
                        quals[i] = rsq;
                    } else if alleles[i] == AlleleType::Ambiguous {
                        // we are already ambiguous, so quality should be 0
                    } else {
                        // check for ambiguity
                        if alleles[i] == rsa {
                            // quals won't always match in local mode
                            // it's also possible to have one quality from global and one from local; lets change this to max
                            quals[i] = quals[i].max(rsq);
                            assert!(quals[i] > 0);
                        } else {
                            // they don't match, change to ambiguous and clear the quality cost
                            alleles[i] = AlleleType::Ambiguous;
                            quals[i] = 0;
                        }
                    }
                }
            }
        }

        // now just send it to the new function
        Self::new(read_name, alleles, quals)
    }

    pub fn read_name(&self) -> &str {
        &self.read_name
    }

    /// Returns the allele value for the given index
    pub fn allele(&self, index: usize) -> AlleleType {
        if self.region.contains(&index) {
            self.alleles[index - self.region.start]
        } else {
            AlleleType::NoOverlap
        }
    }

    /// Returns the quality value for the given index
    pub fn qual(&self, index: usize) -> u8 {
        if self.region.contains(&index) {
            self.quals[index - self.region.start]
        } else {
            0
        }
    }

    /// Returns the range of this segment
    pub fn region(&self) -> &std::ops::Range<usize> {
        &self.region
    }

    /// Returns the number of alleles that are set (i.e. non-ambiguous and overlapping, so 0 or 1)
    pub fn get_num_set(&self) -> usize {
        self.alleles.iter()
            .filter(|&v| *v < AlleleType::Ambiguous)
            .count()
    }

    /// Given a haplotype, this will score the read against that haplotype.
    /// If a haplotype has a 2, no cost is associated with that allele.
    /// # Arguments
    /// `haplotype` - the full haplotype to score, must have the same length as the block
    pub fn score_haplotype(&self, haplotype: &[AlleleType]) -> u64 {
        // with clipped read segments, this no longer holds
        // assert_eq!(self.alleles.len(), haplotype.len());
        // however, the last allele should always be less than haplotype len, or else they've given us something too short
        assert!(self.region.end <= haplotype.len());

        self.score_partial_haplotype(haplotype, 0)
    }

    /// Given a partial haplotype, this will score the read against that haplotype.
    /// The offset values is an index to where to start in our alleles for scoring.
    /// For example, if offset = 10, then alleles[10..] will be compared to haplotype[0..]
    /// If a haplotype has a 2, no cost is associated with that allele.
    /// # Arguments
    /// * `haplotype` - the partial haplotype to score
    /// * `offset` - the offset into the read segment to start scoring
    pub fn score_partial_haplotype(&self, haplotype: &[AlleleType], offset: usize) -> u64 {
        // assertion does not work with clipped segments
        // assert!(haplotype.len()+offset <= self.alleles.len());
        if haplotype.len() + offset <= self.region.start || offset >= self.region.end {
            // the haplotype starts and ends before our first allele, OR
            // the haplotype starts after our last allele, SO
            // return 0 in either case, because there is no overlaps to score
            0
        } else {
            // the minimum comparison is either the first allele OR the offset, whichever is greater
            let min_compare = self.region.start.max(offset);

            // the maximum comparison is either the last allele we have OR the end of the haplotype+start offset
            let max_compare = (self.region.end).min(offset+haplotype.len());
            
            // iterate over the region and add up the quality costs when they mismatch
            (min_compare..max_compare)
                .filter_map(|index| {
                    let allele = self.allele(index);
                    
                    // keep if the allele is non-ambiguous
                    if haplotype[index - offset] < AlleleType::Ambiguous && allele != haplotype[index - offset] {
                        Some(self.qual(index) as u64)
                    } else {
                        None
                    }
                })
                .sum()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constructor() {
        let rs = ReadSegment::new(
            "read_name".to_string(),
            vec![3_u8, 0, 1, 0, 0, 1, 2, 2, 3, 3].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect(),
            vec![0, 1, 2, 3, 4, 5, 6, 7, 0, 0]
        );

        assert_eq!(rs, ReadSegment {
            read_name: "read_name".to_string(),
            alleles: vec![0, 1, 0, 0, 1].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect(),
            quals: vec![1, 2, 3, 4, 5],
            region: 1..6
        });
    }

    #[test]
    fn test_score_haplotype() {
        let rs = ReadSegment::new(
            "read_name".to_string(),
            vec![3_u8, 0, 1, 0, 0, 1, 2, 1, 3, 3].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect(),
            vec![0, 1, 2, 3, 4, 5, 6, 7, 0, 0]
        );
        assert_eq!(rs.region, 1..8);
        assert_eq!(rs.get_num_set(), 6);

        //identical except for missing value in rs
        let haplotype: Vec<AlleleType> = vec![0_u8, 0, 1, 0, 0, 1, 1, 1, 0, 0].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect();
        assert_eq!(rs.score_haplotype(&haplotype), 6);
        
        //fully empty haplotype
        let haplotype: Vec<AlleleType> = vec![2_u8; 10].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect();
        assert_eq!(rs.score_haplotype(&haplotype), 0);

        //fully wrong haplotype
        let haplotype: Vec<AlleleType> = vec![1_u8, 1, 0, 1, 1, 0, 0, 0, 1, 1].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect();
        assert_eq!(rs.score_haplotype(&haplotype), 1+2+3+4+5+6+7);
    }

    #[test]
    fn test_score_partial_haplotype() {
        let rs = ReadSegment::new(
            "read_name".to_string(),
            vec![2_u8, 0, 1, 0, 0, 1, 2, 1, 2, 2].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect(),
            vec![0, 1, 2, 3, 4, 5, 6, 7, 0, 0]
        );

        //identical except for missing value in rs
        let haplotype: Vec<AlleleType> = vec![0, 1, 0, 0, 1, 1, 1].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect();
        assert_eq!(rs.score_partial_haplotype(&haplotype, 1), 6);
        
        //fully empty haplotype
        let haplotype: Vec<AlleleType> = vec![2; 7].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect();
        assert_eq!(rs.score_partial_haplotype(&haplotype, 2), 0);

        //fully wrong haplotype
        let haplotype: Vec<AlleleType> = vec![1, 0, 1, 1, 0, 0, 0].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect();
        assert_eq!(rs.score_partial_haplotype(&haplotype, 1), 1+2+3+4+5+6+7);
        
        for x in 0..haplotype.len() {
            assert_eq!(rs.score_partial_haplotype(&haplotype[x..], 1+x), ((x+1)..8).sum::<usize>() as u64);
        }
    }

    #[test]
    fn test_collapse() {
        let rs1 = ReadSegment::new(
            "read_name".to_string(),
            vec![3, 1, 0, 2, 1, 3, 3].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect(),
            vec![0, 2, 1, 0, 2, 0, 0]
        );
        let rs2 = ReadSegment::new(
            "read_name".to_string(),
            vec![3, 3, 0, 1, 0, 1, 1].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect(),
            vec![0, 0, 1, 2, 2, 1, 1]
        );
        let expected = ReadSegment::new(
            "read_name".to_string(),
            vec![3, 1, 0, 2, 2, 1, 1].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect(),
            vec![0, 2, 1, 0, 0, 1, 1]
        );

        // make sure normal collapsing works
        let collapsed = ReadSegment::collapse(&[rs1.clone(), rs2.clone()]);
        assert_eq!(expected, collapsed);
        assert_eq!(collapsed.region, 1..7);

        // make sure scoring works fine with the 3s present
        //              vec![3, 1, 0, 2, 2, 1, 1]
        let haplotype: Vec<AlleleType> = vec![0_u8, 1, 0, 0, 0, 1, 0].into_iter().map(|v| AlleleType::from_repr(v).unwrap()).collect();
        assert_eq!(collapsed.score_haplotype(&haplotype), 1);

        // check stupid collapsing also
        let collapsed = ReadSegment::collapse(&[rs1.clone()]);
        assert_eq!(collapsed, rs1);
    }
}