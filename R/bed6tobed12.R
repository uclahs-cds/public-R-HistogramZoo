#' Takes a GRanges object and creates a BED12 GRanges object
#'
#' @param gr A GRanges object consisting of the ranges that need to be merged into 1 line of a BED12 file
#'
#' @return A GRanges object of length 1 with mcols matching the additional columns of a BED12 file
#'
#' @examples \dontrun{
#' gr = GenomicRanges::GRanges(
#' seqnames = "chr1",
#' IRanges::IRanges(start = c(1, 4), end = c(3, 6)),
#' strand = "+")
#'
#' bed6tobed12(gr)
#' }
#'
#' @export
bed6tobed12 = function(
  gr
){
  # Generating BED12 columns
  bed12.gr = range(gr)
  S4Vectors::mcols(bed12.gr)$blockCount = length(gr)
  blockSizes = IRanges::width(gr) - 1
  S4Vectors::mcols(bed12.gr)$blockSizes = paste0(blockSizes, ",", collapse = "")
  blockStarts = GenomicRanges::start(gr) - GenomicRanges::start(bed12.gr)
  S4Vectors::mcols(bed12.gr)$blockStarts = paste0(blockStarts, ",", collapse = "")
  bed12.gr
}
