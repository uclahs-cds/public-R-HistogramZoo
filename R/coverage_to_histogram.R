#' Generates a list of Histogram objects from a coverage RLE object and a GRangesList
#'
#' @param regions A GRangesList object of regions to generate Histogram objects
#' @param coverage An RLE object representing genome coverage
#' @param histogram_bin_size An integer representing the width of Histogram bins
coverage_to_histogram = function(
  regions,
  coverage,
  histogram_bin_size
){

  # Initializing
  histogram_coverage = vector("list", length(regions))
  names(histogram_coverage) = names(regions)

  # Calculating coverage
  for(i in seq_along(regions)){
    x = regions[i]
    x_name = names(x)
    x = unlist(x)
    names(x) = NULL
    bins = GenomicRanges::tile(x = x, width = histogram_bin_size)
    bins = unlist(bins)
    GenomeInfoDb::seqlevels(bins) = GenomeInfoDb::seqlevels(coverage)
    cvg = GenomicRanges::binnedAverage(
      bins = bins,
      numvar = coverage,
      varname = "cvg")
    histogram_coverage[[x_name]] <- new_GenomicHistogram(
      histogram_data = cvg$cvg,
      interval_start = GenomicRanges::start(cvg),
      interval_end = GenomicRanges::end(cvg),
      region_id = x_name,
      chr = as.character(GenomicRanges::seqnames(cvg))[1],
      strand = as.character(GenomicRanges::strand(cvg))[1]
    )
  }
  return(histogram_coverage)
}
