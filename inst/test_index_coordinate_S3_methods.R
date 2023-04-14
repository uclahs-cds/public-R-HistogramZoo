

# Preamble ----------------------------------------------------------------
# Writing a set of tests to test S3 methods for fitting histograms
# Updated functions

# Coordinate based functions
# 1. fit_distributions
# 2. fit_uniform
# 3. identify_uniform_segment

# Index based functions
# 4. find_local_optima
# 5. FTC
# 6. meaningful_gaps_local, find_all_meaningful_gap

# Plotting function -------------------------------------------------------


# Fit distributions -------------------------------------------------------

# Numeric
x <- c(1, 2, 3, 4, 5, 4, 3, 2, 1, 1)
  
num_res <- fit_distributions(
  x, 
  dist = "norm", 
  truncated = F, 
  metric = c("mle", "jaccard")
)
  
  # Check that intervals are valid
  
# Histogram (non-1 start)
bins <- index_to_start_end(
  seq(10, 30, 2),
  right = FALSE
)

x <- Histogram(
  x, 
  interval_start = bins$start - 1, 
  interval_end = bins$end
)
  
histogram_res <- fit_distributions(x)
  
  # GenomicHistogram
x <- GenomicHistogram(
  x,
  interval_start = bins$start,
  interval_end = bins$end
)
  
genomic_histogram_res <- fit_distributions(x)
  
