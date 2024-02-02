
use crate::sequence_alignment::edit_distance;

use log::trace;
use std::cmp::Ordering;

/// All the variant types we are currently allowing
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq, PartialOrd, Ord)]
pub enum VariantType {
    /// REF and ALT are both length = 1
    Snv=0,
    /// REF length = 1, ALT length > 1
    Insertion,
    /// REF length > 1, ALT length = 1
    Deletion,
    /// REF and ALT lengths > 1
    Indel,
    /// Must have two alleles and be tagged with SVTYPE=INS; ALT >= REF
    SvInsertion,
    /// Must have two alleles and be tagged with SVTYPE=DEL; ALT <= REF
    SvDeletion,
    /// Must have two alleles and be tagged with SVTYPE=DUP
    SvDuplication,
    /// Must have two alleles and be tagged with SVTYPE=INV
    SvInversion,
    /// Must have two alleles and be tagged with SVTYPE=BND
    SvBreakend,
    /// Must have two alleles and be tagged with TRID=####
    TandemRepeat,
    /// Something that doesn't match the above criteria, must be 1 or 2 alleles
    Unknown // make sure Unknown is always the last one in the list
}

/// Zygosity definitions, mostly used elsewhere
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq, PartialOrd, Ord)]
pub enum Zygosity {
    HomozygousReference=0,
    Heterozygous,
    HomozygousAlternate,
    Unknown // make sure Unknown is always the last one in the list
}

#[derive(thiserror::Error, Debug)]
pub enum VariantError {
    #[error("allele0 length must match reference length when index_allele0=0")]
    Allele0RefLen,
    #[error("allele{index} must be length 1")]
    AlleleLen1{ index: usize },
    #[error("allele{index} is empty (length = 0)")]
    EmptyAllele{ index: usize },
    #[error("index_allele0 must be < index_allele1")]
    IndexAlleleOrder,
    #[error("{variant_type:?} does not support multi-allelic sites")]
    MultiAllelicNotAllowed { variant_type: VariantType },
    #[error("reference must have length > 1")]
    RefLenGT1,
    #[error("SV deletion ALT length must be <= REF length")]
    SvDeletionLen,
    #[error("SV insertion ALT length must be >= REF length")]
    SvInsertionLen
}

/// A variant definition structure.
/// It currently assumes that chromosome is fixed and that the variant is a SNP.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Variant {
    /// The vcf index from the input datasets
    vcf_index: usize,
    /// The type of variant represented by this entry
    variant_type: VariantType,
    /// The coordinate of the event in the VCF file, 0-based
    position: i64,
    /// The length of the reference allele
    ref_len: usize,
    /// the start position of the allele (0-based), will be <= position
    prefix_len: usize,
    /// the end position (0-based, exclusive), will be >= position+ref_len
    postfix_len: usize,
    /// the first allele value
    allele0: Vec<u8>,
    /// the second allele value
    allele1: Vec<u8>,

    //these only matter for multi-allelic sites; usize is "proper" type, but u8 will be nice and compact
    /// the index of allele0, typically 0 (REF)
    index_allele0: u8,
    /// the index of allele1, typically 1 (usually len(ALT) == 1, so it's 1)
    index_allele1: u8,

    // auxiliary booleans
    /// if true, flags this is a variant to ignore for _some_ reason
    is_ignored: bool
}

impl Variant {
    /// Creates a new single-nucleotide variant (SNV).
    /// For SNV variants, all alleles must be exactly 1 bp long.
    /// # Arguments
    /// * `vcf_index` - the index of the source VCF file
    /// * `position` - the coordinate of the variant in a contig
    /// * `allele0` - the first allele (usually REF)
    /// * `allele1` - the second allele (usually ALT[0])
    /// * `index_allele0` - the index for allele0, typically 0 for REF
    /// * `index_allele1` - the index for allele1, typically 1 for simple heterozygous variants
    /// # Panics
    /// * if `index_allele0 > index_allele1`
    /// * if the provided sequences do not match a single-nucleotide variant
    pub fn new_snv(vcf_index: usize, position: i64, allele0: Vec<u8>, allele1: Vec<u8>, index_allele0: u8, index_allele1: u8) -> Result<Variant, VariantError> {
        // we always assume alleles come "sorted" and they are heterozygous
        if index_allele0 >= index_allele1 {
            return Err(VariantError::IndexAlleleOrder);
        }

        // SNV alleles must be length 1
        if allele0.len() != 1 {
            return Err(VariantError::AlleleLen1 { index: 0 });
        }
        if allele1.len() != 1 {
            return Err(VariantError::AlleleLen1 { index: 1 });
        }

        Ok(Variant {
            vcf_index,
            variant_type: VariantType::Snv,
            position,
            ref_len: 1,
            prefix_len: 0,
            postfix_len: 0,
            allele0,
            allele1,
            index_allele0,
            index_allele1,
            is_ignored: false
        })
    }

    /// Creates a new deletion variant.
    /// Deletions must have a REF allele long than 1 bp, and all ALT alleles must be exactly 1 bp long.
    /// # Arguments
    /// * `vcf_index` - the index of the source VCF file
    /// * `position` - the coordinate of the variant in a contig
    /// * `ref_len` - the length of the reference allele
    /// * `allele0` - the first allele (usually REF)
    /// * `allele1` - the second allele (usually ALT[0])
    /// * `index_allele0` - the index for allele0, typically 0 for REF
    /// * `index_allele1` - the index for allele1, typically 1 for simple heterozygous variants
    /// # Panics
    /// * if `index_allele0 > index_allele1`
    /// * if the provided sequences do not match a deletion variant
    /// * if the reference allele is passed in and it does not have the same length as `ref_len`
    pub fn new_deletion(vcf_index: usize, position: i64, ref_len: usize, allele0: Vec<u8>, allele1: Vec<u8>, index_allele0: u8, index_allele1: u8) -> Result<Variant, VariantError> {
        // we always assume alleles come "sorted" and they are heterozygous
        if index_allele0 >= index_allele1 {
            return Err(VariantError::IndexAlleleOrder);
        }
        
        // reference length must be greater than 1 to be a deletion
        if ref_len <= 1 {
            return Err(VariantError::RefLenGT1);
        }
        
        if index_allele0 == 0 {
            // this allele is also the reference allele
            if allele0.len() != ref_len {
                return Err(VariantError::Allele0RefLen);
            }
        } else {
            // this allele is not the reference, must be a multi-allelic site; but all deletion alts have len = 1
            if allele0.len() != 1 {
                return Err(VariantError::AlleleLen1 { index: 0 });
            }
        }

        // this one must always be length 1
        if allele1.len() != 1{
            return Err(VariantError::AlleleLen1 { index: 1 });
        }
        
        // make sure the alleles start with the same thing
        if allele0[0] != allele1[0] {
            /*
            Counter example to requiring alleleX[0] be equal; this is rare, but it does seem to happen
            chr12	117794450	.	ACACACCAACATGCACACT	T
            */
            trace!("Deletion alleles are unexpected: {position}, {ref_len}, {allele0:?}, {allele1:?}");
        }
        Ok(Variant {
            vcf_index,
            variant_type: VariantType::Deletion,
            position,
            ref_len,
            prefix_len: 0,
            postfix_len: 0,
            allele0,
            allele1,
            index_allele0,
            index_allele1,
            is_ignored: false
        })
    }

    /// Creates a new insertion variant.
    /// Insertions must have a REF allele exactly 1 bp long, and all ALT alleles must be longer than 1 bp.
    /// # Arguments
    /// * `vcf_index` - the index of the source VCF file
    /// * `position` - the coordinate of the variant in a contig
    /// * `allele0` - the first allele (usually REF)
    /// * `allele1` - the second allele (usually ALT[0])
    /// * `index_allele0` - the index for allele0, typically 0 for REF
    /// * `index_allele1` - the index for allele1, typically 1 for simple heterozygous variants
    /// # Panics
    /// * if `index_allele0 > index_allele1`
    /// * if the provided sequences do not match an insertion variant
    pub fn new_insertion(vcf_index: usize, position: i64, allele0: Vec<u8>, allele1: Vec<u8>, index_allele0: u8, index_allele1: u8) -> Result<Variant, VariantError> {
        // we always assume alleles come "sorted" and they are heterozygous
        if index_allele0 >= index_allele1 {
            return Err(VariantError::IndexAlleleOrder);
        }
        
        if index_allele0 == 0 {
            // if reference allele is present, it must be length 1 for this type
            if allele0.len() != 1 {
                return Err(VariantError::AlleleLen1 { index: 0 });
            }
        } else {
            // allele0 isn't reference, so it must be >= 1 due to multi-allelics
            //     chr1	2122634	.	T	C,TG	14.1
            if allele0.is_empty() {
                return Err(VariantError::EmptyAllele { index: 0 });
            }
        }
        // we have to do >= because of some multi-allelics:
        //     chr1	286158	.	A	ATG,G	34.4
        if allele1.is_empty() {
            return Err(VariantError::EmptyAllele { index: 1 });
        }

        // make sure the alleles start with the same thing
        if allele0[0] != allele1[0] {
            // no counter example searched for yet, but probably exists, we'll leave this trace for now
            trace!("Insertion alleles are unexpected: {position}, 1, {allele0:?}, {allele1:?}");
        }
        Ok(Variant {
            vcf_index,
            variant_type: VariantType::Insertion,
            position,
            ref_len: 1,
            prefix_len: 0,
            postfix_len: 0,
            allele0,
            allele1,
            index_allele0,
            index_allele1,
            is_ignored: false
        })
    }

    /// Creates a new indel variant.
    /// All indels alleles must be more than 1 bp long.
    /// # Arguments
    /// * `vcf_index` - the index of the source VCF file
    /// * `position` - the coordinate of the variant in a contig
    /// * `ref_len` - the length of the reference allele
    /// * `allele0` - the first allele (usually REF)
    /// * `allele1` - the second allele (usually ALT[0])
    /// * `index_allele0` - the index for allele0, typically 0 for REF
    /// * `index_allele1` - the index for allele1, typically 1 for simple heterozygous variants
    /// # Panics
    /// * if `index_allele0 > index_allele1`
    /// * if the provided sequences do not match an indel variant
    /// * if the reference allele is passed in and it does not have the same length as `ref_len`
    pub fn new_indel(vcf_index: usize, position: i64, ref_len: usize, allele0: Vec<u8>, allele1: Vec<u8>, index_allele0: u8, index_allele1: u8) -> Result<Variant, VariantError> {
        // we always assume alleles come "sorted" and they are heterozygous
        if index_allele0 >= index_allele1 {
            return Err(VariantError::IndexAlleleOrder);
        }

        // reference length must be greater than 1 to be an indel, but ALTs can really be any length after that (>=1 anyways)
        if ref_len <= 1 {
            return Err(VariantError::RefLenGT1);
        }
        
        if index_allele0 == 0 {
            // this allele is also the reference allele
            if allele0.len() != ref_len {
                return Err(VariantError::Allele0RefLen);
            }
        } else {
            // it's not a reference allele, since this is an indel, length can be anything >= 1
            if allele0.is_empty() {
                return Err(VariantError::EmptyAllele { index: 0 });
            }
        }

        // this one just has to be >= 1
        if allele1.is_empty() {
            return Err(VariantError::EmptyAllele { index: 1 });
        }
        
        // there's no real reason to believe in any shared sequence between alleles
        // we've seen it not work above, not worth even trying to codify warning here IMO
        // assert!(???)
        
        Ok(Variant {
            vcf_index,
            variant_type: VariantType::Indel,
            position,
            ref_len,
            prefix_len: 0,
            postfix_len: 0,
            allele0,
            allele1,
            index_allele0,
            index_allele1,
            is_ignored: false
        })
    }

    /// Creates a new SV deletion variant.
    /// SV deletions must have a REF allele long than 1 bp, and all ALT alleles must be exactly 1 bp long.
    /// # Arguments
    /// * `vcf_index` - the index of the source VCF file
    /// * `position` - the coordinate of the variant in a contig
    /// * `ref_len` - the length of the reference allele
    /// * `allele0` - the first allele (usually REF)
    /// * `allele1` - the second allele (usually ALT[0])
    /// * `index_allele0` - the index for allele0, typically 0 for REF
    /// * `index_allele1` - the index for allele1, typically 1 for simple heterozygous variants
    /// # Panics
    /// * if `index_allele0 > index_allele1`
    /// * if the provided sequences do not match a deletion variant
    /// * if the reference allele is passed in and it does not have the same length as `ref_len`
    pub fn new_sv_deletion(vcf_index: usize, position: i64, ref_len: usize, allele0: Vec<u8>, allele1: Vec<u8>, index_allele0: u8, index_allele1: u8) -> Result<Variant, VariantError> {
        // we always assume alleles come "sorted" and they are heterozygous
        if index_allele0 >= index_allele1 {
            return Err(VariantError::IndexAlleleOrder);
        }

        // this is one difference from plain Deletion, must be 0 and 1 in VCF
        if index_allele0 != 0 || index_allele1 != 1 {
            return Err(VariantError::MultiAllelicNotAllowed { variant_type: VariantType::SvDeletion })
        }
        
        // this allele is also the reference allele
        if allele0.len() != ref_len {
            return Err(VariantError::Allele0RefLen);
        }

        // restriction lifted such that now allele0 must be >= allele1
        if allele1.len() > allele0.len() {
            return Err(VariantError::SvDeletionLen);
        }
        
        // allele1 is always <= allele0, so just make sure it's not empty
        if allele1.is_empty() {
            return Err(VariantError::EmptyAllele { index: 1 });
        }
        
        // make sure the alleles start with the same thing
        if allele0[0] != allele1[0] {
            /*
            Counter example to requiring alleleX[0] be equal; this is rare, but it does seem to happen
            chr12	117794450	.	ACACACCAACATGCACACT	T
            */
            trace!("Deletion alleles are unexpected: {}, {}, {:?}, {:?}", position, ref_len, allele0, allele1);
        }
        Ok(Variant {
            vcf_index,
            variant_type: VariantType::SvDeletion,
            position,
            ref_len,
            prefix_len: 0,
            postfix_len: 0,
            allele0,
            allele1,
            index_allele0,
            index_allele1,
            is_ignored: false
        })
    }

    /// Creates a new SV insertion variant.
    /// SV insertions must have a REF allele exactly 1 bp long, and all ALT alleles must be longer than 1 bp.
    /// # Arguments
    /// * `vcf_index` - the index of the source VCF file
    /// * `position` - the coordinate of the variant in a contig
    /// * `ref_len` - the length of the reference allele
    /// * `allele0` - the first allele (usually REF)
    /// * `allele1` - the second allele (usually ALT[0])
    /// * `index_allele0` - the index for allele0, typically 0 for REF
    /// * `index_allele1` - the index for allele1, typically 1 for simple heterozygous variants
    /// # Panics
    /// * if `index_allele0 > index_allele1`
    /// * if the provided sequences do not match an insertion variant
    pub fn new_sv_insertion(vcf_index: usize, position: i64, ref_len: usize, allele0: Vec<u8>, allele1: Vec<u8>, index_allele0: u8, index_allele1: u8) -> Result<Variant, VariantError> {
        // we always assume alleles come "sorted" and they are heterozygous
        if index_allele0 >= index_allele1 {
            return Err(VariantError::IndexAlleleOrder);
        }
        
        // this is one difference from plain Insertion, must be 0 and 1 in VCF
        if index_allele0 != 0 || index_allele1 != 1 {
            return Err(VariantError::MultiAllelicNotAllowed { variant_type: VariantType::SvInsertion })
        }
        
        // this allele is also the reference allele
        if allele0.len() != ref_len {
            return Err(VariantError::Allele0RefLen);
        }

        // restriction lifted such that now allele1 must be >= allele0
        if allele1.len() < allele0.len() {
            return Err(VariantError::SvInsertionLen);
        }
        
        // allele0 is always <= allele1, so just make sure it's not empty
        if allele0.is_empty() {
            return Err(VariantError::EmptyAllele { index: 0 });
        }

        // make sure the alleles start with the same thing
        if allele0[0] != allele1[0] {
            // no counter example searched for yet, but probably exists, we'll leave this trace for now
            trace!("Insertion alleles are unexpected: {}, {}, {:?}, {:?}", position, 1, allele0, allele1);
        }
        Ok(Variant {
            vcf_index,
            variant_type: VariantType::SvInsertion,
            position,
            ref_len,
            prefix_len: 0,
            postfix_len: 0,
            allele0,
            allele1,
            index_allele0,
            index_allele1,
            is_ignored: false
        })
    }

    /// Creates a new tandem repeat variant, functionally these act very similar to indel types.
    /// All tandem repeat alleles must be at least 1 bp long by VCF definition.
    /// # Arguments
    /// * `vcf_index` - the index of the source VCF file
    /// * `position` - the coordinate of the variant in a contig
    /// * `ref_len` - the length of the reference allele
    /// * `allele0` - the first allele (usually REF)
    /// * `allele1` - the second allele (usually ALT[0])
    /// * `index_allele0` - the index for allele0, typically 0 for REF
    /// * `index_allele1` - the index for allele1, typically 1 for simple heterozygous variants
    /// # Panics
    /// * if `index_allele0 > index_allele1`
    /// * if the provided sequences do not match a tandem repeat variant
    /// * if the reference allele is passed in and it does not have the same length as `ref_len`
    pub fn new_tandem_repeat(vcf_index: usize, position: i64, ref_len: usize, allele0: Vec<u8>, allele1: Vec<u8>, index_allele0: u8, index_allele1: u8) -> Result<Variant, VariantError> {
        // we always assume alleles come "sorted" and they are heterozygous
        if index_allele0 >= index_allele1 {
            return Err(VariantError::IndexAlleleOrder);
        }

        // all alleles must be >= 1 for tandem repeats, most are longer though
        if allele0.is_empty() {
            return Err(VariantError::EmptyAllele { index: 0 });
        }
        if allele1.is_empty() {
            return Err(VariantError::EmptyAllele { index: 1 });
        }
        
        if index_allele0 == 0 {
            // this allele is also the reference allele, make sure the lengths match
            if allele0.len() != ref_len {
                return Err(VariantError::Allele0RefLen);
            }
        } else {
            // it's not a reference allele, so nothing to check here
        }
        
        Ok(Variant {
            vcf_index,
            variant_type: VariantType::TandemRepeat,
            position,
            ref_len,
            prefix_len: 0,
            postfix_len: 0,
            allele0,
            allele1,
            index_allele0,
            index_allele1,
            is_ignored: false
        })
    }

    /// This will add a prefix to each allele, generally reference genome sequence that will allow for better matching.
    /// # Arguments
    /// * `prefix` - the sequence to pre-pend to each allele
    pub fn add_reference_prefix(&mut self, prefix: &[u8]) {
        // make sure we don't set our reference start coordinate to less than 0
        let prefix_len: usize = prefix.len();
        assert!(prefix_len <= self.position as usize - self.prefix_len);
        
        // allele0, pre-pend is basically copy
        let mut new_allele0: Vec<u8> = Vec::with_capacity(self.allele0.len()+prefix_len);
        new_allele0.extend_from_slice(prefix);
        new_allele0.extend_from_slice(&self.allele0);
        self.allele0 = new_allele0;

        // same for allele1
        let mut new_allele1: Vec<u8> = Vec::with_capacity(self.allele1.len()+prefix_len);
        new_allele1.extend_from_slice(prefix);
        new_allele1.extend_from_slice(&self.allele1);
        self.allele1 = new_allele1;

        // finally, adjust the start coordinates
        self.prefix_len += prefix_len;
    }

    /// This will add a postfix to each allele, generally reference genome sequence that will allow for better matching.
    /// # Arguments
    /// * `postfix` - the sequence to append to each allele
    pub fn add_reference_postfix(&mut self, postfix: &[u8]) {
        // easier operation, just extend the existing vecs
        self.allele0.extend_from_slice(postfix);
        self.allele1.extend_from_slice(postfix);

        // finally, adjust the end coordinates
        self.postfix_len += postfix.len();
    }

    /// This will trim the postfix down to a smaller size.
    pub fn truncate_reference_postfix(&mut self, truncate_amount: usize) {
        // sanity check that we are only truncating the postfix
        assert!(truncate_amount <= self.postfix_len);

        // truncate the alleles and shrink the postfix size
        self.allele0.truncate(self.allele0.len() - truncate_amount);
        self.allele1.truncate(self.allele1.len() - truncate_amount);
        self.postfix_len -= truncate_amount;
    }

    pub fn get_vcf_index(&self) -> usize {
        self.vcf_index
    }

    pub fn get_type(&self) -> VariantType {
        self.variant_type
    }

    pub fn position(&self) -> i64 {
        self.position
    }

    pub fn get_ref_len(&self) -> usize {
        self.ref_len
    }

    pub fn get_prefix_len(&self) -> usize {
        self.prefix_len
    }

    pub fn get_postfix_len(&self) -> usize {
        self.postfix_len
    }

    pub fn get_allele0(&self) -> &[u8] {
        &self.allele0
    }

    pub fn get_allele1(&self) -> &[u8] {
        &self.allele1
    }

    pub fn is_ignored(&self) -> bool {
        self.is_ignored
    }

    pub fn set_ignored(&mut self) {
        self.is_ignored = true;
    }

    pub fn get_truncated_allele0(&self) -> &[u8] {
        let start: usize = self.prefix_len;
        let end: usize = self.allele0.len() - self.postfix_len;
        &self.allele0[start..end]
    }

    pub fn get_truncated_allele1(&self) -> &[u8] {
        let start: usize = self.prefix_len;
        let end: usize = self.allele1.len() - self.postfix_len;
        &self.allele1[start..end]
    }

    /// This will determine the best matching allele (0 or 1) or return 2 if neither match.
    /// Primary purpose of this is to convert all variant observations into a 0/1 scheme.
    /// This method requires an exact match of the allele.
    /// # Arguments
    /// * `allele` - the allele that needs to get converted to a 0 or 1 (or 2 if neither match)
    pub fn match_allele(&self, allele: &[u8]) -> u8 {
        if allele == &self.allele0[..] {
            0
        } else if allele == &self.allele1[..] {
            1
        } else {
            2
        }
    }

    /// This will determine the closest matching allele (0 or 1) based on edit distance, or return 2 if they are equi-distant.
    /// This method does not require an exact match of the alleles.
    /// Returns a tuple of the (allele chosen, min edit distance, other edit distance).
    /// # Arguments
    /// * `allele` - the allele sequence to compare to our internal alleles
    pub fn closest_allele(&self, allele: &[u8]) -> (u8, usize, usize) {
        self.closest_allele_clip(allele, 0, 0)
    }

    /// This will determine the closest matching allele (0 or 1) based on edit distance, or return 2 if they are equi-distant.
    /// This method does not require an exact match of the alleles, and allows for you to clip bases on the internal allele sequence.
    /// This is most useful when you have to clip the provided allele to to incomplete matching.
    /// Returns a tuple of the (allele chosen, min edit distance, other edit distance).
    /// # Arguments
    /// * `allele` - the allele sequence to compare to our internal alleles
    /// * `offset` - will skip this many bases internally for calculating edit distance
    pub fn closest_allele_clip(&self, allele: &[u8], head_clip: usize, tail_clip: usize) -> (u8, usize, usize) {
        assert!(head_clip <= self.prefix_len);
        assert!(tail_clip <= self.postfix_len);
        let d0: usize = edit_distance(allele, &self.allele0[head_clip..(self.allele0.len() - tail_clip)]);
        let d1: usize = edit_distance(allele, &self.allele1[head_clip..(self.allele1.len() - tail_clip)]);
        trace!("clipping: {} {}", head_clip, tail_clip);
        trace!("obs{:?}", allele);
        trace!("a0 {:?} => {}", &self.allele0[head_clip..(self.allele0.len() - tail_clip)], d0);
        trace!("a1 {:?} => {}", &self.allele1[head_clip..(self.allele1.len() - tail_clip)], d1);
        match d0.cmp(&d1) {
            // d0 is less, return that
            Ordering::Less => (0, d0, d1),
            // d1 is less, return that
            Ordering::Greater => (1, d1, d0),
            // equidistant, so undetermined
            Ordering::Equal => (2, d0, d1)
        }
    }

    /// This will return the index allele for a given haplotype index.
    /// Input must always be 0 or 1, but it might get converted to something else at multi-allelic sites.
    /// # Arguments
    /// * `index` - must be 0, 1, or 2 (unknown)
    /// # Panics
    /// * if anything other than 0, 1, or 2 is provided
    pub fn convert_index(&self, index: u8) -> u8 {
        if index == 0 {
            self.index_allele0
        } else if index == 1 {
            self.index_allele1
        } else if index == 2 {
            // we just need some indicator that it's undetermined, this will work for now
            u8::MAX  
        } else {
            panic!("index must be 0, 1, or 2");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_basic_snv() {
        let variant = Variant::new_snv(
            0, 1,
            b"A".to_vec(), b"C".to_vec(),
            0, 1
        ).unwrap();
        assert_eq!(variant.get_type(), VariantType::Snv);
        assert_eq!(variant.position(), 1);
        assert_eq!(variant.get_ref_len(), 1);
        assert_eq!(variant.match_allele(b"A"), 0);
        assert_eq!(variant.match_allele(b"C"), 1);
        assert_eq!(variant.match_allele(b"G"), 2);
        assert_eq!(variant.match_allele(b"T"), 2);
        assert_eq!(variant.convert_index(0), 0);
        assert_eq!(variant.convert_index(1), 1);
        assert_eq!(variant.convert_index(2), u8::MAX);
    }

    #[test]
    fn test_basic_deletion() {
        // this is the deletion we mostly expect
        let variant = Variant::new_deletion(
            0, 10, 3,
            b"AGT".to_vec(), b"A".to_vec(),
            0, 1
        ).unwrap();
        assert_eq!(variant.get_type(), VariantType::Deletion);
        assert_eq!(variant.position(), 10);
        assert_eq!(variant.get_ref_len(), 3);
        assert_eq!(variant.match_allele(b"AGT"), 0);
        assert_eq!(variant.match_allele(b"A"), 1);
        assert_eq!(variant.match_allele(b"AG"), 2);

        // multi-allelic deletion, must still be length 1 though
        let variant = Variant::new_deletion(
            0, 10, 4,
            b"C".to_vec(), b"A".to_vec(),
            1, 2
        ).unwrap();
        assert_eq!(variant.get_type(), VariantType::Deletion);
        assert_eq!(variant.position(), 10);
        assert_eq!(variant.get_ref_len(), 4);
        assert_eq!(variant.match_allele(b"ACCC"), 2);
        assert_eq!(variant.match_allele(b"C"), 0);
        assert_eq!(variant.match_allele(b"A"), 1);
        assert_eq!(variant.convert_index(0), 1);
        assert_eq!(variant.convert_index(1), 2);
    }

    #[test]
    fn test_basic_insertion() {
        let variant = Variant::new_insertion(
            0, 20,
            b"A".to_vec(), b"AGT".to_vec(),
            0, 1
        ).unwrap();
        assert_eq!(variant.get_type(), VariantType::Insertion);
        assert_eq!(variant.position(), 20);
        assert_eq!(variant.get_ref_len(), 1);
        assert_eq!(variant.match_allele(b"A"), 0);
        assert_eq!(variant.match_allele(b"AGT"), 1);
        assert_eq!(variant.match_allele(b"AG"), 2);
    }

    #[test]
    fn test_basic_indel() {
        // models AG -> A / AGT
        let variant = Variant::new_indel(
            0, 20, 2,
            b"A".to_vec(), b"AGT".to_vec(),
            1, 2
        ).unwrap();
        assert_eq!(variant.get_type(), VariantType::Indel);
        assert_eq!(variant.position(), 20);
        assert_eq!(variant.get_ref_len(), 2);
        assert_eq!(variant.match_allele(b"A"), 0);
        assert_eq!(variant.match_allele(b"AGT"), 1);
        assert_eq!(variant.match_allele(b"AG"), 2);
    }

    #[test]
    fn test_sv_insertion() {
        let variant = Variant::new_sv_insertion(
            0, 20, 1,
            b"A".to_vec(), b"AGT".to_vec(),
            0, 1
        ).unwrap();
        assert_eq!(variant.get_type(), VariantType::SvInsertion);
        assert_eq!(variant.position(), 20);
        assert_eq!(variant.get_ref_len(), 1);

        // TODO: replace this with the matching we will do with SVs
        assert_eq!(variant.match_allele(b"A"), 0);
        assert_eq!(variant.match_allele(b"AGT"), 1);
        assert_eq!(variant.match_allele(b"AG"), 2);
    }

    #[test]
    fn test_sv_deletion() {
        let variant = Variant::new_sv_deletion(
            0, 10, 3,
            b"AGT".to_vec(), b"A".to_vec(),
            0, 1
        ).unwrap();
        assert_eq!(variant.get_type(), VariantType::SvDeletion);
        assert_eq!(variant.position(), 10);
        assert_eq!(variant.get_ref_len(), 3);

        // TODO: replace this with the matching we will do with SVs
        assert_eq!(variant.match_allele(b"AGT"), 0);
        assert_eq!(variant.match_allele(b"A"), 1);
        assert_eq!(variant.match_allele(b"AG"), 2);
    }

    #[test]
    fn test_tandem_repeat() {
        let variant = Variant::new_tandem_repeat(
            0, 10, 4, 
            b"AAAC".to_vec(),
            b"AAACAAAC".to_vec(),
            0, 1
        ).unwrap();
        assert_eq!(variant.get_type(), VariantType::TandemRepeat);
        assert_eq!(variant.position(), 10);
        assert_eq!(variant.get_ref_len(), 4);

        assert_eq!(variant.match_allele(b"AAAC"), 0);
        assert_eq!(variant.match_allele(b"AAACAAAC"), 1);
        assert_eq!(variant.match_allele(b"AAACAA"), 2);
    }

    #[test]
    fn test_reference_adjustment() {
        // models AG -> A / AGT
        let mut variant = Variant::new_indel(
            0, 20, 2,
            b"A".to_vec(), b"AGT".to_vec(),
            1, 2
        ).unwrap();

        // make sure no fixins yet
        assert_eq!(variant.get_prefix_len(), 0);
        assert_eq!(variant.get_postfix_len(), 0);
        
        let prefix: Vec<u8> = b"AC".to_vec();
        variant.add_reference_prefix(&prefix);
        let postfix: Vec<u8> = b"GGCC".to_vec();
        variant.add_reference_postfix(&postfix);

        assert_eq!(variant.get_truncated_allele0(), b"A");
        assert_eq!(variant.get_truncated_allele1(), b"AGT");

        // trims off the extra 'C' we added
        variant.truncate_reference_postfix(1);
        
        // make sure nothing here changes
        assert_eq!(variant.get_type(), VariantType::Indel);
        assert_eq!(variant.position(), 20);
        assert_eq!(variant.get_ref_len(), 2);

        // check this new stuff
        assert_eq!(variant.get_prefix_len(), 2);
        assert_eq!(variant.get_postfix_len(), 3);

        // original alleles will not match exactly anymore
        assert_eq!(variant.match_allele(b"A"), 2);
        assert_eq!(variant.match_allele(b"AGT"), 2);
        assert_eq!(variant.match_allele(b"AG"), 2);

        // inexact without the reference data will return weird results
        assert_eq!(variant.closest_allele(b"A"), (0, 5, 7));
        assert_eq!(variant.closest_allele(b"AGT"), (0, 4, 5));
        assert_eq!(variant.closest_allele(b"AG"), (0, 4, 6));

        // now lets inexact with the extensions
        assert_eq!(variant.closest_allele(b"ACAGGC"), (0, 0, 2));
        assert_eq!(variant.closest_allele(b"ACAGTGGC"), (1, 0, 2));
        assert_eq!(variant.closest_allele(b"ACAGGGC"), (2, 1, 1));
    }
}
