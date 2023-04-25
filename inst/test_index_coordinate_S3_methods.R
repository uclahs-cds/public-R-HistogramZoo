

library(HistogramZoo)
library(BoutrosLab.plotting.general)

# Preamble ----------------------------------------------------------------
# Writing a set of tests to test S3 methods for fitting histograms
# Updated functions

# Issue with shift on numeric data

# Coordinate based functions
# 1. fit_distributions
# - The issue with gamma is obvious - we need an additional offset
# - Unclear issue with bin_width
# 2. fit_uniform
# 3. identify_uniform_segment

# Index based functions
# 4. find_local_optima
# 5. FTC
# 6. meaningful_gaps_local, find_all_meaningful_gap

# Plotting function -------------------------------------------------------

distribution_colours <- c(
  "norm" = "darkorange",
  "gamma" = "chartreuse4",
  "gamma_flip" = "chartreuse3",
  "unif" = "darkorchid4"
)

metric_lty <- c(
  "mle" = 1,
  "jaccard" = 2,
  "intersection" = 3,
  "chisq" = 4,
  "ks" = 5,
  "mse" = 6
)

plot_fit_distributions <- function(res, x, empirical_data){
  
  plotting_data <- do.call(rbind.data.frame, lapply(res, function(res_i){
    data.frame(
      "x" = x,
      "y" = res_i$dens(x, mpar = res_i$par, scale = T),
      "metric" = res_i[['metric']],
      "distribution" = res_i[['dist']],
      "tag" = paste0(res_i[['metric']], "-", res_i[['dist']])
    )
  }))
  
  
  plt <- create.scatterplot(
    y ~ x,
    plotting_data,
    # plotting data
    groups = plotting_data$tag,
    type = "l",
    col = distribution_colours[unique(plotting_data$distribution)],
    lty = metric_lty[unique(plotting_data$metric)],
    # Axes
    xlab.cex = 1,
    xlab.label = "x",
    ylab.cex = 1,
    ylab.label = "Data",
    xaxis.cex = 1,
    yaxis.cex = 1,
    xlimits = c(min(x), max(x)),
    # Adding points
    add.points = T,
    points.x = x,
    points.y = empirical_data,
    points.col = "red"
  )
  
  return(plt)
}

find_midpoint <- function(start, end){
  (start + end)/2
}

# Fit distributions -------------------------------------------------------

# Numeric
x <- seq(1, 10)
y <- c(1, 2, 3, 4, 5, 4, 3, 2, 1, 1)
  
num_res <- fit_distributions(y, metric = "jaccard", dist = "norm")

plot_fit_distributions(
  res = num_res, 
  x = x,
  empirical_data = y
)

# Histogram (non-1 start)
bins <- index_to_start_end(
  # seq(10, 20, 1),
  seq(10, 30, 2),
  right = FALSE
)

histogram_non_zero_start <- Histogram(
  y, 
  interval_start = bins$start - 1, 
  interval_end = bins$end
)
  
histogram_res <- fit_distributions(histogram_non_zero_start) # , metric = "jaccard", dist = "norm")

plot_fit_distributions(
  res = histogram_res[1:5], 
  x = find_midpoint(histogram_non_zero_start$interval_start, histogram_non_zero_start$interval_end),
  empirical_data = y
)

# GenomicHistogram
genomic_histogram_non_zero_start <- GenomicHistogram(
  y,
  interval_start = bins$start,
  interval_end = bins$end
)
  
genomic_histogram_res <- fit_distributions(genomic_histogram_non_zero_start, dist = "norm") # , metric = "jaccard", dist = "norm")

plot_fit_distributions(
  res = genomic_histogram_res, 
  x = find_midpoint(genomic_histogram_non_zero_start$consecutive_start, genomic_histogram_non_zero_start$consecutive_end),
  empirical_data = y
)


# Fit uniform -------------------------------------------------------------

# Numeric
x <- seq(1, 10)
y_unif <- rep(1, 10) + runif(10)

num_res <- fit_uniform(y_unif)

plot_fit_distributions(
  res = list(num_res), 
  x = x,
  empirical_data = y_unif
)

# Histogram
unif_bins <- index_to_start_end(
  # seq(10, 20, 1),
  seq(10, 30, 2),
  right = FALSE
)

histogram_non_zero_start <- Histogram(
  y_unif, 
  interval_start = unif_bins$start - 1, 
  interval_end = unif_bins$end
)

histogram_res <- fit_uniform(histogram_non_zero_start) 

plot_fit_distributions(
  res = list(histogram_res), 
  x = find_midpoint(histogram_non_zero_start$interval_start, histogram_non_zero_start$interval_end),
  empirical_data = y_unif
)

# GenomicHistogram
genomic_histogram_non_zero_start <- GenomicHistogram(
  y_unif,
  interval_start = unif_bins$start,
  interval_end = unif_bins$end
)

genomic_histogram_res <- fit_uniform(genomic_histogram_non_zero_start)

plot_fit_distributions(
  res = list(genomic_histogram_res), 
  x = find_midpoint(genomic_histogram_non_zero_start$consecutive_start, genomic_histogram_non_zero_start$consecutive_end),
  empirical_data = y_unif
)

# Identify uniform segment ------------------------------------------------

# Numeric
x <- seq(1, 10)
y_unif_trim <- c(1, 1, rep(2, 7), 1)

num_res <- identify_uniform_segment(y_unif_trim)

plot_fit_distributions(
  res = list(num_res), 
  x = x,
  empirical_data = y_unif_trim
)

# Histogram
unif_bins <- index_to_start_end(
  # seq(10, 20, 1),
  seq(10, 30, 2),
  right = FALSE
)

histogram_non_zero_start <- Histogram(
  y_unif_trim, 
  interval_start = unif_bins$start - 1, 
  interval_end = unif_bins$end
)

histogram_res <- identify_uniform_segment(histogram_non_zero_start) 

plot_fit_distributions(
  res = list(histogram_res), 
  x = find_midpoint(histogram_non_zero_start$interval_start, histogram_non_zero_start$interval_end),
  empirical_data = y_unif_trim
)

# GenomicHistogram
genomic_histogram_non_zero_start <- GenomicHistogram(
  y_unif_trim,
  interval_start = unif_bins$start,
  interval_end = unif_bins$end
)

genomic_histogram_res <- identify_uniform_segment(genomic_histogram_non_zero_start)

plot_fit_distributions(
  res = list(genomic_histogram_res), 
  x = find_midpoint(genomic_histogram_non_zero_start$consecutive_start, genomic_histogram_non_zero_start$consecutive_end),
  empirical_data = y_unif_trim
)
