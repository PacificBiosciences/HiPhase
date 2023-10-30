
// This was all provided by Daniel Baker as a way to get VCF indexing since rust_htslib does not have a user-friendly method yet.
// Parts have been tweaked for readability.

/// The error type we can generate from trying to index
#[derive(Debug)]
pub struct BcfBuildError {
    pub msg: String,
}

impl BcfBuildError {
    pub fn error_message(error: i32) -> &'static str {
        match error {
           -1 => "indexing failed",
           -2 => "opening @fn failed",
           -3 => "format not indexable",
           -4 => "failed to create and/or save the index",
            _ => "unknown error",
        }
    }
}
impl std::error::Error for BcfBuildError {}

impl std::fmt::Display for BcfBuildError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "BcfBuildError{{msg: {}}}", self.msg)
    }
}

/// Build a bcf or vcf.gz index.
/// Builds tbi or csi depending on if build_tbi is set.
pub fn build_bcf_index<P: AsRef<std::path::Path>>(
    bcf_path: P,
    idx_path: Option<P>,
    n_threads: u32,
    build_tbi: bool,
) -> Result<(), BcfBuildError> {
    let min_shift = if build_tbi {0} else {14};
    let idx_path_cstr = idx_path.map(|p| rust_htslib::utils::path_to_cstring(&p).expect("path_to_cstring"));
    let ret = unsafe {
        rust_htslib::htslib::bcf_index_build3(
            rust_htslib::utils::path_to_cstring(&bcf_path).unwrap().as_ptr(),
            idx_path_cstr.map_or(std::ptr::null(), |p| p.as_ptr()),
            min_shift,
            n_threads as i32,
        )
    };
    match ret {
        0 => Ok(()),
        e => Err(BcfBuildError {
            msg: format!("Failed to build  bcf index. Error: {e:?}/{}", BcfBuildError::error_message(e)),
        }),
    }
}
