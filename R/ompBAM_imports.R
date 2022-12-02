#' @useDynLib ompBench, .registration = TRUE
#' @import Rcpp
#' @import zlibbioc
NULL

#' @export
idxstats_ompBAM <- function(bam_file, n_threads, verbose = TRUE) {
    idxstats_pbam(bam_file, n_threads, verbose)
}

#' @export
idxstats_htslib <- function(bam_file, n_threads, verbose = TRUE) {
    idxstats_hts(bam_file, n_threads, verbose)
}

#' @export
idxstats_htslib_omp <- function(bam_file, n_threads, verbose = TRUE) {
    idxstats_hts_omp(bam_file, n_threads, verbose)
}