#' Finds the largest uniform segment that is longer than threshold
#'
#' @param x numeric vector representing the density of a histogram
#' @param metric one of `jaccard`, `intersection`, `ks`, `mse`, `chisq` indicating metrics to use for fit optimization
#' @param threshold numeric, indicating the minimum proportion of the subsegment which should be tested
#' @param stepsize integer, indicating the stepsize (relative to the histogram bins) to take in the search for the uniform subsegment
#' @param max_sd_size numeric, the number of standard deviations of the computed metric distribution away from the optimal uniform which
#' has maximum length
#'
#' @return a list representing a uniform model with the following data
#' \describe{
#'     \item{par}{A character string denoting the region_id of the Histogram}
#'     \item{dist}{The distribution name}
#'     \item{metric}{The metric used to fit the distribution}
#'     \item{value}{The fitted value of the metric function}
#'     \item{dens}{A function that returns the density of the fitted distribution}
#'     \item{seg_start}{start index of the interval}
#'     \item{seg_end}{end index of the interval}
#' }
#'
#' @rdname identify_uniform_segment
#' @export
identify_uniform_segment <- function(
    x,
    metric = c("jaccard", "intersection", "ks", "mse", "chisq"),
    threshold = 0.5,
    stepsize = 1,
    max_sd_size = 1
){
  UseMethod("identify_uniform_segment")
}

#' @rdname identify_uniform_segment
#' @exportS3Method identify_uniform_segment numeric
identify_uniform_segment.numeric <- function(
    x,
    metric = c("jaccard", "intersection", "ks", "mse", "chisq"),
    threshold = 0.5,
    stepsize = 1,
    max_sd_size = 1
){
  interval_midpoint <- seq(1, length(x), 1)
  identify_uniform_segment_helper(
    x = x,
    interval_start = interval_midpoint - 0.5,
    interval_end = interval_midpoint + 0.5,
    interval_midpoint = interval_midpoint,
    metric = metric,
    threshold = threshold,
    stepsize = stepsize,
    max_sd_size = max_sd_size
  )
}

#' @rdname identify_uniform_segment
#' @exportS3Method identify_uniform_segment Histogram
identify_uniform_segment.Histogram <- function(
    x,
    metric = c("jaccard", "intersection", "ks", "mse", "chisq"),
    threshold = 0.5,
    stepsize = 1,
    max_sd_size = 1
){
  identify_uniform_segment_helper(
    x = x$histogram_data,
    interval_start = x$interval_start,
    interval_end = x$interval_end,
    interval_midpoint = find_midpoint(x),
    metric = metric,
    threshold = threshold,
    stepsize = stepsize,
    max_sd_size = max_sd_size
  )
}

#' @rdname identify_uniform_segment
#' @exportS3Method identify_uniform_segment GenomicHistogram
identify_uniform_segment.GenomicHistogram <- function(
    x,
    metric = c("jaccard", "intersection", "ks", "mse", "chisq"),
    threshold = 0.5,
    stepsize = 1,
    max_sd_size = 1
){
  identify_uniform_segment_helper(
    x = x$histogram_data,
    interval_start = x$consecutive_start - 0.5,
    interval_end = x$consecutive_end + 0.5,
    interval_midpoint = find_midpoint(x),
    metric = metric,
    threshold = threshold,
    stepsize = stepsize,
    max_sd_size = max_sd_size
  )
}

identify_uniform_segment_helper <- function(
  x,
  interval_start,
  interval_end,
  interval_midpoint,
  metric = c("jaccard", "intersection", "ks", "mse", "chisq"),
  threshold = 0.5,
  stepsize = 1,
  max_sd_size = 1
){

  # Error checking
  metric <- match.arg(metric)
  if(!is.numeric(threshold) | length(threshold) != 1 ){
    stop("threshold must be a numeric of length 1")
  }
  if(!(threshold >= 0 & threshold <= 1)) {
    stop("threshold must be between 0 and 1")
  }
  if(!is_equal_integer(stepsize) | !(stepsize > 0) | length(stepsize) != 1){
    stop("stepsize must be functional as a positive integer of length 1")
  }
  if(!is.numeric(max_sd_size) | !(max_sd_size >= 0) | length(max_sd_size) != 1){
    stop("max_sd_size must be a positive or zero numeric of length 1")
  }

  # Set-up
  num_bins <- length(x)
  min_seg_size <- ceiling(num_bins * threshold)

  res <- lapply(seq(from = 1, to = num_bins - min_seg_size, by = stepsize), function(a) {
    lapply(seq(from = min_seg_size + a, to = num_bins, by = stepsize), function(b) {
        fitted_unif <- fit_uniform_helper(
            x = x[a:b],
            interval_start = interval_start[a:b],
            interval_end = interval_end[a:b],
            interval_midpoint = interval_midpoint[a:b],
            metric = metric
          )
        fitted_unif$seg_start <- a
        fitted_unif$seg_end <- b
        return(fitted_unif)
    })
  })

  # Extracting stats
  res <- unlist(res, recursive = F)
  res_df <- do.call(rbind.data.frame, lapply(res, `[`, c('seg_start', 'seg_end', 'value')))
  res_df$index <- 1:nrow(res_df)

  # Calculating length
  res_df$length <- res_df$seg_end - res_df$seg_start + 1

  # Correct for Jaccard/Intersection
  res_df$value <- if (metric %in% c("jaccard", "intersection")) (1 - res_df$value) else res_df$value

  # Select the longest interval that is within XX sd of the maximum
  min_metric <- min(res_df$value, na.rm = T)
  sd_metric <- stats::sd(res_df$value, na.rm = T)
  sd_metric <- ifelse(is.na(sd_metric), 0, sd_metric) # If there's only 1 case, set it to 0
  res_sd_range <- res_df[res_df$value <= min_metric + sd_metric * max_sd_size, ]
  max_interval_index <- res_sd_range$index[which.max(res_sd_range$length)]

  return(
    res[[max_interval_index]]
  )

}
