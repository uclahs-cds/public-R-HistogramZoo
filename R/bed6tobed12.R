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
  S4Vectors::mcols(bed12.gr)$score = 0
  S4Vectors::mcols(bed12.gr)$thickStart = GenomicRanges::start(bed12.gr)
  S4Vectors::mcols(bed12.gr)$thickEnd = GenomicRanges::end(bed12.gr)
  S4Vectors::mcols(bed12.gr)$itemRgb = 0
  S4Vectors::mcols(bed12.gr)$blockCount = length(gr)
  blockSizes = IRanges::width(gr) - 1
  S4Vectors::mcols(bed12.gr)$blockSizes = paste0(blockSizes, ",", collapse = "")
  blockStarts = GenomicRanges::start(gr) - GenomicRanges::start(bed12.gr)
  S4Vectors::mcols(bed12.gr)$blockStarts = paste0(blockStarts, ",", collapse = "")
  bed12.gr
}
