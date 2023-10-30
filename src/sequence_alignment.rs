
/// Returns the edit distance between two u8 Vecs by doing the full grid calculation.
/// This version is row-based (rows are length of v1) for the main loop.
/// # Arguments
/// * `v1` - the first sequence
/// * `v2` - the second sequence
pub fn edit_distance(v1: &[u8], v2: &[u8]) -> usize {
    // structured such that each "row" is the length of v1 (i.e. v1 is conceptually on the x-axis)
    let l1: usize = v1.len();
    let mut row: Vec<usize> = vec![0; l1+1];
    let mut prev_row: Vec<usize> = (0..l1+1).collect();
    
    // go through each row
    for (i, &c2) in v2.iter().enumerate() {
        row[0] = i+1;
        for (j, &c1) in v1.iter().enumerate() {
            row[j+1] = [
                // skip a character in v2
                prev_row[j+1]+1,
                // skip a character in v1 
                row[j]+1,
                // diagonal match/mismatch
                prev_row[j]+({
                    if c1 == c2 {
                        0
                    } else {
                        1
                    }
                })
            ].into_iter().min().unwrap();
        }

        // swap the rows at the end of each iteration
        std::mem::swap(&mut row, &mut prev_row);
    }

    prev_row[l1]
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_edit_distance() {
        let v1: Vec<u8> = vec![0, 1, 2, 4, 5];
        let v2: Vec<u8> = vec![0, 1, 3, 4, 5];
        let v3: Vec<u8> = vec![1, 2, 3, 5];
        let v4: Vec<u8> = vec![];

        assert_eq!(edit_distance(&v1, &v1), 0);
        assert_eq!(edit_distance(&v1, &v2), 1);
        assert_eq!(edit_distance(&v1, &v3), 2);
        assert_eq!(edit_distance(&v1, &v4), 5);

        assert_eq!(edit_distance(&v2, &v2), 0);
        assert_eq!(edit_distance(&v2, &v3), 3);
        assert_eq!(edit_distance(&v2, &v4), 5);
        
        assert_eq!(edit_distance(&v3, &v3), 0);
        assert_eq!(edit_distance(&v3, &v4), 4);

        assert_eq!(edit_distance(&v4, &v4), 0);
    }

    #[test]
    fn test_edit_error_001() {
        let v1 = [65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 67, 65, 65, 65];
        let v2 = [65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 67, 65, 65, 65, 65, 65, 65, 67, 65, 65, 65];
        let v3 = [65, 65, 65, 65, 65, 65, 65, 65, 65, 65,     65, 65, 65, 65, 65, 65, 67, 65, 65, 65];

        assert_eq!(edit_distance(&v1, &v3), 1);
        assert_eq!(edit_distance(&v2, &v3), 1);
        assert_eq!(edit_distance(&v3, &v1), 1);
        assert_eq!(edit_distance(&v3, &v2), 1);
    }
}