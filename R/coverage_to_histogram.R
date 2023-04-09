
#' A helper function for coverage_to_histogram
#'
#' @param split_coverage A GRangesList object where each element is a GRanges object representing 1 bin
#' and each element in the GRanges object contains an mcols column `cvg` representing the coverage for
#' that element
#'
#' return A GRanges object where each element is one bin and `cvg` represents the weighted average of the
#' coverage of its components
weighted_coverage <- function(split_coverage){
  seg_length <- IRanges::width(split_coverage)
  gr <- base::range(split_coverage)
  gr$cvg <- sum( split_coverage$cvg * seg_length) / sum(seg_length)
  return(gr)
}

#' Generates a GenomicHistogram object from coverage data
#'
#' @param coverage An RLE object representing genome coverage
#' @param region a GRanges object representing the region to be binned for 1 GenomicHistogram
#' @param region_id character, region_id of the resulting histogram
#' @param histogram_bin_size An integer representing the width of Histogram bins
#'
#' @return a GenomicHistogram object
#'
#' @export
coverage_to_histogram = function(
    region,
    region_id,
    coverage,
    histogram_bin_size
){

  region_range <- base::range(region)
  introns <- GenomicRanges::setdiff(region_range, region)
  bins <- unlist(GenomicRanges::tile(x = region, width = 1))
  breaks <- seq(1, length(bins), histogram_bin_size)
  breaks <- unique(c(breaks, length(bins)))

  if(histogram_bin_size == 1){
    break_start <- breaks
    break_end <- breaks
  } else if(length(breaks) == 2){
    break_start <- 1
    break_end <- breaks[2]
  } else {
    break_start <- breaks[1:(length(breaks)-1)]
    break_end <- c(breaks[2:(length(breaks)-1)]-1, tail(breaks, n = 1))
  }

  bins <- lapply(seq_along(break_start), function(i) IRanges::reduce(bins[break_start[i]:break_end[i]]))
  bins <- GenomicRanges::GRangesList(bins)
  names(bins) <- generate_grangeslist_identifiers(bins)
  bins <- unlist(bins)
  GenomeInfoDb::seqlevels(bins) <- GenomeInfoDb::seqlevels(coverage)

  cvg <- GenomicRanges::binnedAverage(
    bins = bins,
    numvar = coverage,
    varname = "cvg"
  )

  cvg <- GenomicRanges::split(cvg, f = names(cvg))
  cvg[sapply(cvg, length) > 1] <- GenomicRanges::GRangesList(sapply(cvg[sapply(cvg, length) > 1], weighted_coverage))
  cvg <- unlist(cvg)
  cvg <- sort(cvg)

  # Generating new histogram
  return(
    HistogramZoo:::new_GenomicHistogram(
      histogram_data = unname(cvg$cvg),
      interval_start = GenomicRanges::start(cvg),
      interval_end = GenomicRanges::end(cvg),
      bin_width = as.integer(histogram_bin_size),
      region_id = region_id,
      chr = as.character(GenomicRanges::seqnames(cvg))[1],
      strand = as.character(GenomicRanges::strand(cvg))[1],
      intron_start = GenomicRanges::start(introns),
      intron_end = GenomicRanges::end(introns)
    )
  )
}
