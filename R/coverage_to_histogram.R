#' Generates a GenomicHistogram object from coverage data
#'
#' @param coverage An RLE object representing genome coverage
#' @param region a GRanges object representing the region to be binned for 1 GenomicHistogram
#' @param region_id character, region_id of the resulting histogram
#' @param histogram_bin_size An integer representing the width of Histogram bins
#'
#' @return a GenomicHistogram object
coverage_to_histogram = function(
  region,
  region_id,
  coverage,
  histogram_bin_size
){

  # Generating bins
  bins <- unlist(
    GenomicRanges::tile(x = region, width = histogram_bin_size)
  )
  GenomeInfoDb::seqlevels(bins) <- GenomeInfoDb::seqlevels(coverage)

  # Computing coverage
  cvg <- GenomicRanges::binnedAverage(
    bins = bins,
    numvar = coverage,
    varname = "cvg"
  )

  # Generating new histogram
  return(
    new_GenomicHistogram(
      histogram_data = cvg$cvg,
      interval_start = GenomicRanges::start(cvg),
      interval_end = GenomicRanges::end(cvg),
      bin_width = as.integer(histogram_bin_size),
      region_id = region_id,
      chr = as.character(GenomicRanges::seqnames(cvg))[1],
      strand = as.character(GenomicRanges::strand(cvg))[1]
    )
  )
}
