
#' A helper function for coverage_to_histogram that returns the weighted coverage of a disjoint bin
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


#' A helper function for coverage_to_histogram that returns a dataframe from a set of ranges
#'
#' @param parent A GRange object from which to extract chromosome and strand
#' @param ranges_dataframe A dataframe with `start` and `end` columns representing genomic coordinates for the
#' new GRanges object
#'
#' return A GRanges object with the chromosome and strand of the parent and the ranges of the ranges_dataframe
generate_granges_from_parent_and_ranges <- function(parent, ranges_dataframe){
  GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(parent)[1],
    IRanges::IRanges(
      start = ranges_dataframe$start,
      end = ranges_dataframe$end
    ),
    strand = GenomicRanges::strand(parent)[1]
  )
}

#' A helper function for coverage_to_histogram that returns a dataframe from a set of ranges
#'
#' @param coverage An RLE object representing genome coverage
#' @param bins A GRanges object representing a set of bins
#'
#' return A GRanges object a "cvg" column representing computed coverage
compute_coverage_on_bins <- function(coverage, bins){
  GenomeInfoDb::seqlevels(bins) <- GenomeInfoDb::seqlevels(coverage)
  GenomicRanges::binnedAverage(bins = bins, numvar = coverage, varname = "cvg")
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
coverage_to_histogram <- function(
    region,
    region_id,
    coverage,
    histogram_bin_size
){

  # Initializing
  region_range <- base::range(region)
  introns <- GenomicRanges::setdiff(region_range, region)

  if(histogram_bin_size == 1){
    bins <- unlist(GenomicRanges::tile(x = region, width = 1))
    cvg <- compute_coverage_on_bins(coverage, bins)

  } else {
  
    # Identify start and end coordinates
    region_start <- GenomicRanges::start(region)
    region_end <- GenomicRanges::end(region)
    index <- lapply(seq_along(region_start), function(i) seq(region_start[i], region_end[i], 1))
    index <- do.call(c, index)
    breaks <- seq(1, length(index), histogram_bin_size)
    breaks <- unique(c(breaks, length(index)))
    breaks <- index[breaks]

    # Generate a table
    df <- index_to_start_end(breaks, right = FALSE)

    # If introns do not exist
    if(length(introns) == 0){
      bins <- generate_granges_from_parent_and_ranges(parent = region, ranges_dataframe = df)
      cvg <- compute_coverage_on_bins(coverage, bins)
    } else { # Aiyah, introns and split bins!

      intron_start <- GenomicRanges::start(introns)
      intron_end <- GenomicRanges::end(introns)
      split_bin_info <- apply(df, 1, function(x) {
        idx <- intron_start >= x['start'] & intron_end <= x['end']
        list(any(idx), intron_start[idx], intron_end[idx])
        }
      )
      split_bin_idx <- sapply(split_bin_info, `[[`, 1)
      split_bins <- df[split_bin_idx,]
      split_bin_start <- sapply(split_bin_info, `[[`, 2)[split_bin_idx]
      split_bin_end <- sapply(split_bin_info, `[[`, 3)[split_bin_idx]

      split_bin_coverage <- lapply(seq_along(split_bin_start), function(i){
        # This isn't a perfectly efficient use of remove_max_gaps but we trust there aren't *that* many introns
        bins_wo_introns <- remove_max_gaps(split_bins[i,], data.frame("start" = split_bin_start[[i]], "end" = split_bin_end[[i]]))
        bins_gr <- generate_granges_from_parent_and_ranges(region, bins_wo_introns)
        bins_coverage <- compute_coverage_on_bins( coverage, bins_gr)
        weighted_coverage(bins_coverage)
      })

      continuous_bins <- df[!split_bin_idx,]
      continuous_gr <- generate_granges_from_parent_and_ranges(region, continuous_bins)
      continuous_coverage <- compute_coverage_on_bins( coverage, continuous_gr)

      cvg <- do.call(c, c(split_bin_coverage, list(continuous_coverage)))
      cvg <- sort(cvg)
    }
  }

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
