
library(HistogramZoo)
library(GenomicRanges)

# Preamble ----------------------------------------------------------------
# Refactoring coverage to histogram to address the following objectives
# Last bin can be shorter (like Histogram)


# Function ----------------------------------------------------------------

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
    HistogramZoo:::new_GenomicHistogram(
      histogram_data = cvg$cvg,
      interval_start = GenomicRanges::start(cvg),
      interval_end = GenomicRanges::end(cvg),
      bin_width = as.integer(histogram_bin_size),
      region_id = region_id,
      chr = as.character(GenomicRanges::seqnames(cvg))[1],
      strand = as.character(GenomicRanges::strand(cvg))[1],
      # The following two lines are the only changes to the original function
      # as a result of the addition of introns
      intron_start = integer(),
      intron_end = integer()
    )
  )
}

generate_grangeslist_identifiers <- function(gr_list){
  gr_ranges <- range(gr_list)
  return(
    paste0(
      GenomicRanges::seqnames(gr_ranges), ":",
      GenomicRanges::start(gr_ranges), "-",
      GenomicRanges::end(gr_ranges), ":",
      GenomicRanges::strand(gr_ranges)
    )
  )
}

weighted_coverage <- function(split_coverage){
  seg_length <- IRanges::width(split_coverage)
  gr <- base::range(split_coverage)
  gr$cvg <- sum( split_coverage$cvg * seg_length) / sum(seg_length)
  return(gr)
}

# Test Data ---------------------------------------------------------------

one_segment_region = GenomicRanges::GRanges(
  seqnames = "chr1",
  IRanges::IRanges(
    start = c(1),
    end = c(10)
  ),
  strand = "+"
)

two_segment_region = GenomicRanges::GRanges(
  seqnames = "chr1",
  IRanges::IRanges(
    start = c(1, 21),
    end = c(10, 30)
  ),
  strand = "+"
)

multi_intron_region = GenomicRanges::GRanges(
  seqnames = "chr1",
  IRanges::IRanges(
    start = c(1, 6, 9),
    end = c(4, 7, 20)
  ),
  strand = "+"
)

# coverage
coverage <- GenomicRanges::coverage(
  GenomicRanges::GRanges(
    seqnames = "chr1",
    IRanges::IRanges(
      start = c(1, 21, 31),
      end = c(20, 30, 40)
    ),
    strand = "+"
  )
)

# uneven coverage
uneven_coverage <- GenomicRanges::coverage(
  GenomicRanges::GRanges(
    seqnames = "chr1",
    IRanges::IRanges(
      start = c(1, 21, 21),
      end = c(20, 30, 30)
    ),
    strand = "+"
  )
)

# region id
region_id = "TEST"

# Fixing tiling -----------------------------------------------------------

coverage_to_histogram_new = function(
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

# benchmark ---------------------------------------------------------------

compare_functions <- function(region, coverage, region_id, histogram_bin_size){
  
  start_time <- Sys.time()
  old = coverage_to_histogram(
    region = region,
    coverage = coverage,
    region_id =  region_id,
    histogram_bin_size = histogram_bin_size
  )
  end_time <- Sys.time()
  t1 = end_time - start_time
  
  start_time <- Sys.time()
  new = coverage_to_histogram_new(
      region = region,
      coverage = coverage,
      region_id =  region_id,
      histogram_bin_size = histogram_bin_size
  )
  end_time <- Sys.time()
  t2 = end_time - start_time
  
  cat("Slow down: ", as.numeric(t2)/as.numeric(t1), "\n")
  
  # return(list("New" = t2, "Old" = t1))
  return(new)
}


# Test Cases --------------------------------------------------------------

# Bin size 1 - base case
base_case <- compare_functions(
  region = one_segment_region,
  coverage = coverage,
  region_id =  region_id,
  histogram_bin_size = 1
)

# Bin size 2 - base case
base_case_two <- compare_functions(
  region = one_segment_region,
  coverage = coverage,
  region_id =  region_id,
  histogram_bin_size = 2
)

# Checking that the last bin is indeed a different size if length of object is different
non_divisible_bin_length <- compare_functions(
  region = one_segment_region,
  coverage = coverage,
  region_id =  region_id,
  histogram_bin_size = 4
)

# Checking that bins span introns
base_intron_case <- compare_functions(
  region = two_segment_region,
  coverage = coverage,
  region_id =  region_id,
  histogram_bin_size = 4
)

# Multiple introns per bin
multi_intron_case <- compare_functions(
  region = multi_intron_region,
  coverage = coverage,
  region_id =  region_id,
  histogram_bin_size = 4
)

# Uneven coverage
uneven_coverage_case <- compare_functions(
  region = two_segment_region,
  coverage = uneven_coverage,
  region_id =  region_id,
  histogram_bin_size = 4
)