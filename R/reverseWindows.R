#' Reverse the effective orientation of selected chromosomes
#'
#' @description
#' When comparing two assemblies, homologous chromosomes can be assembled in
#' opposite orientations. As a result, the per-window coverage and GC vectors of
#' such chromosomes run in opposite directions along the genome.
#'
#' `reverseWindows()` flips the effective orientation of the requested
#' chromosomes so that windows can be compared position by position. It keeps
#' the window indexes (`iid`, `start`, `stop`, and the chromosome identifier)
#' unchanged, and instead reverses the value vectors (coverage, GC-corrected
#' coverage when present, and the per-window validity flags) within each
#' selected chromosome. The (modified) object is returned.
#'
#' For an `ISA` object, segmentation is **recomputed** for every
#' donor-recipient pair that contains a reversed chromosome, because these
#' pairs' difference profiles have changed orientation. Segments of unaffected
#' pairs are preserved unchanged.
#'
#' @param x An object of class `RCA` or `ISA`.
#' @param chr_ids Character vector of chromosome IDs in the object's meta table
#'   (`chr_id` column) whose orientation should be reversed. The short display
#'   name (`chr_name`), used mainly for plotting labels, must not be used here.
#' @param ... Further arguments passed to methods.
#'
#' @return The modified object, of the same class as `x`.
#'
#' @examples
#' \dontrun{
#' rev_rca <- reverseWindows(rca, chr_ids = c("OY997261.1", "CP155614.1"))
#' rev_isa <- reverseWindows(isa, chr_ids = "OY997261.1")
#' }
#' @export
reverseWindows <- function(x, chr_ids, ...) {
    UseMethod("reverseWindows")
}