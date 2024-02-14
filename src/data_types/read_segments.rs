
/// Container for a read segment that has been converted into a variant representation
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReadSegment {
    /// the read name
    read_name: String,
    /// the actual alleles for the read, should always be 0, 1, 2 (ambiguous), or 3 (non-overlapping/undefined).
    /// anything other than 0 or 1 are basically ignored in all functions
    alleles: Vec<u8>,
    /// the associated quality values for converting 0 <--> 1; undefined alleles should have qual == 0
    quals: Vec<u8>,
    /// the index of the first defined 0/1 allele, inclusive
    first_allele: usize,
    /// the index of the last defined 0/1 allele, inclusive
    last_allele: usize
}

impl ReadSegment {
    /// Creates a new read segment from a set of alleles and quality values.
    /// Note that the segment should have the same length as the phase block, even if the segment does not actually span the full block.
    /// # Arguments
    /// * `alleles` - the alleles, should all be 0, 1, 2 (ambiguous) or 3 (non-overlapping/undefined); the first and last defined values determine the length of the segment
    /// * `quals` - cost to convert an allele from 0 <--> 1, any undefined/ambiguous alleles should be qual = 0
    /// # Panics
    /// * if `allele.len() != quals.len()`
    pub fn new(read_name: String, alleles: Vec<u8>, quals: Vec<u8>) -> ReadSegment {
        assert_eq!(alleles.len(), quals.len());
        let (first_allele, _) = alleles.iter().enumerate()
            .find(|(_i, &a)| a < 2).unwrap_or((alleles.len(), &2));
        let (last_allele, _) = alleles.iter().enumerate().rev()
            .find(|(_i, &a)| a < 2).unwrap_or((alleles.len(), &2));
        ReadSegment {
            read_name,
            alleles,
            quals,
            first_allele,
            last_allele
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

        let num_alleles: usize = read_segments[0].get_num_alleles();
        let read_name: String = read_segments[0].read_name().to_string();
        let mut alleles: Vec<u8> = vec![3; num_alleles];
        let mut quals: Vec<u8> = vec![0; num_alleles];
        for rs in read_segments.iter() {
            let rs_alleles = rs.alleles();
            let rs_quals = rs.quals();
            assert_eq!(num_alleles, rs.get_num_alleles());
            assert_eq!(read_name, rs.read_name());

            for (i, &rsa) in rs_alleles.iter().enumerate() {
                // if rsa is unset, we skip everything
                if rsa != 3 {
                    if alleles[i] == 3 {
                        alleles[i] = rsa;
                        quals[i] = rs_quals[i];
                    } else if alleles[i] == 2 {
                        // we are already ambiguous, so quality should be 0
                    } else {
                        // check for ambiguity
                        if alleles[i] == rsa {
                            // quals won't always match in local mode
                            // it's also possible to have one quality from global and one from local; lets change this to max
                            quals[i] = quals[i].max(rs_quals[i]);
                            assert!(quals[i] > 0);
                        } else {
                            // they don't match, change to ambiguous
                            alleles[i] = 2;
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

    pub fn get_num_alleles(&self) -> usize {
        self.alleles.len()
    }

    pub fn alleles(&self) -> &[u8] {
        &self.alleles[..]
    }

    pub fn quals(&self) -> &[u8] {
        &self.quals[..]
    }

    pub fn first_allele(&self) -> usize {
        self.first_allele
    }

    pub fn last_allele(&self) -> usize {
        self.last_allele
    }

    /// Returns the range of this segment, e.g. [first_allele..last_allele+1)
    pub fn get_range(&self) -> std::ops::Range<usize> {
        self.first_allele..(self.last_allele+1)
    }

    /// Returns the number of alleles that are set (i.e. non-ambiguous and overlapping, so 0 or 1)
    pub fn get_num_set(&self) -> usize {
        self.alleles.iter()
            .filter(|&v| *v < 2)
            .count()
    }

    /// Given a haplotype, this will score the read against that haplotype.
    /// If a haplotype has a 2, no cost is associated with that allele.
    /// # Arguments
    /// `haplotype` - the full haplotype to score, must have the same length as the block
    pub fn score_haplotype(&self, haplotype: &[u8]) -> u64 {
        assert_eq!(self.alleles.len(), haplotype.len());
        self.score_partial_haplotype(haplotype, 0)
    }

    /// Given a partial haplotype, this will score the read against that haplotype.
    /// The offset values is an index to where to start in our alleles for scoring.
    /// For example, if offset = 10, then alleles[10..] will be compared to haplotype[0..]
    /// If a haplotype has a 2, no cost is associated with that allele.
    /// # Arguments
    /// * `haplotype` - the partial haplotype to score
    /// * `offset` - the offset into the read segment to start scoring
    pub fn score_partial_haplotype(&self, haplotype: &[u8], offset: usize) -> u64 {
        //info!("rs {}+{} <= {}?", haplotype.len(), offset, self.alleles.len());
        assert!(haplotype.len()+offset <= self.alleles.len());
        if haplotype.len() + offset <= self.first_allele || offset > self.last_allele {
            // the haplotype starts and ends before our first allele, OR
            // the haplotype starts after our last allele, SO
            // return 0 in either case, because there is no overlaps to score
            0
        } else {
            // the minimum comparison is either the first allele OR the offset, whichever is greater
            let min_compare = self.first_allele.max(offset);

            // if the allele component is greater, then we need to shift into the haplotype
            let offset_shift = min_compare - offset;

            // the maximum comparison is either the last allele we have OR the end of the haplotype+start offset
            let max_compare = (self.last_allele+1).min(offset+haplotype.len());
            
            self.quals[min_compare..max_compare].iter().enumerate()
                .filter(|(i, _q)| haplotype[*i+offset_shift] < 2 && self.alleles[*i+min_compare] != haplotype[*i+offset_shift])
                .map(|(_i, &q)| q as u64)
                .sum()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_score_haplotype() {
        let rs = ReadSegment::new(
            "read_name".to_string(),
            vec![2, 0, 1, 0, 0, 1, 2, 1, 2, 2],
            vec![0, 1, 1, 1, 1, 1, 1, 1, 0, 0]
        );
        assert_eq!(rs.first_allele, 1);
        assert_eq!(rs.last_allele, 7);
        assert_eq!(rs.get_num_set(), 6);

        //identical except for missing value in rs
        let haplotype = vec![0, 0, 1, 0, 0, 1, 1, 1, 0, 0];
        assert_eq!(rs.score_haplotype(&haplotype), 1);
        
        //fully empty haplotype
        let haplotype = vec![2; 10];
        assert_eq!(rs.score_haplotype(&haplotype), 0);

        //fully wrong haplotype
        let haplotype = vec![1, 1, 0, 1, 1, 0, 0, 0, 1, 1];
        assert_eq!(rs.score_haplotype(&haplotype), 7);
    }

    #[test]
    fn test_score_partial_haplotype() {
        let rs = ReadSegment::new(
            "read_name".to_string(),
            vec![2, 0, 1, 0, 0, 1, 2, 1, 2, 2],
            vec![0, 1, 1, 1, 1, 1, 1, 1, 0, 0]
        );

        //identical except for missing value in rs
        let haplotype = vec![0, 1, 0, 0, 1, 1, 1];
        assert_eq!(rs.score_partial_haplotype(&haplotype, 1), 1);
        
        //fully empty haplotype
        let haplotype = vec![2; 7];
        assert_eq!(rs.score_partial_haplotype(&haplotype, 2), 0);

        //fully wrong haplotype
        let haplotype = vec![1, 0, 1, 1, 0, 0, 0];
        assert_eq!(rs.score_partial_haplotype(&haplotype, 1), 7);
        for x in 0..haplotype.len() {
            assert_eq!(rs.score_partial_haplotype(&haplotype[x..], 1+x), 7-x as u64);
        }
    }

    #[test]
    fn test_collapse() {
        let rs1 = ReadSegment::new(
            "read_name".to_string(),
            vec![3, 1, 0, 2, 1, 3, 3],
            vec![0, 2, 1, 0, 2, 0, 0]
        );
        let rs2 = ReadSegment::new(
            "read_name".to_string(),
            vec![3, 3, 0, 1, 0, 1, 1],
            vec![0, 0, 1, 2, 2, 1, 1]
        );
        let expected = ReadSegment::new(
            "read_name".to_string(),
            vec![3, 1, 0, 2, 2, 1, 1],
            vec![0, 2, 1, 0, 0, 1, 1]
        );

        // make sure normal collapsing works
        let collapsed = ReadSegment::collapse(&[rs1.clone(), rs2.clone()]);
        assert_eq!(expected, collapsed);
        assert_eq!(collapsed.first_allele, 1);
        assert_eq!(collapsed.last_allele, 6);

        // make sure scoring works fine with the 3s present
        //              vec![3, 1, 0, 2, 2, 1, 1]
        let haplotype = vec![0, 1, 0, 0, 0, 1, 0];
        assert_eq!(collapsed.score_haplotype(&haplotype), 1);

        // check stupid collapsing also
        let collapsed = ReadSegment::collapse(&[rs1.clone()]);
        assert_eq!(collapsed, rs1);
    }
}