#' Generates a list of histogram count vectors from a coverage RLE object and a GRangesList
#'
#' @param regions TODO
#' @param coverage.rle TODO
#' @param histogram.bin.size TODO
#'
#' @return TODO
#' @export
coverage.to.histogram = function(
  regions,
  coverage.rle,
  histogram.bin.size
){

  # Initializing
  histogram.coverage = vector("list", length(regions))
  names(histogram.coverage) = names(regions)

  # Calculating coverage
  for(i in seq_along(regions)){
    x = regions[i]
    x.name = names(x)
    x = unlist(x)
    names(x) = NULL
    bins = GenomicRanges::tile(x = x, width = histogram.bin.size)
    bins = unlist(bins)
    GenomeInfoDb::seqlevels(bins) = GenomeInfoDb::seqlevels(coverage.rle)
    cvg = GenomicRanges::binnedAverage(
      bins = bins,
      numvar = coverage.rle,
      varname = "cvg")
    histogram.coverage[[x.name]] <- cvg$cvg
  }
  histogram.coverage
}
