# TODO: Check that this is compatible with Base 0 and Base 1 formatting
#' Takes a GRanges object and creates a BED12 GRanges object
#'
#' @param gr TODO
#'
#' @return
#' @export
bed6tobed12 = function(
  gr
){
  bed12.gr = range(gr)
  # blockCount
  S4Vectors::mcols(bed12.gr)$blockCount = length(gr)
  # blockSizes
  blockSizes = IRanges::width(gr) - 1
  S4Vectors::mcols(bed12.gr)$blockSizes = paste0(blockSizes, ",", collapse = "")
  # blockStarts
  blockStarts = GenomicRanges::start(gr) - GenomicRanges::start(bed12.gr)
  S4Vectors::mcols(bed12.gr)$blockStarts = paste0(blockStarts, ",", collapse = "")
  # Return GRanges object
  bed12.gr
}
