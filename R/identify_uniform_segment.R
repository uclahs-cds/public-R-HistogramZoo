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
#' @export
identify_uniform_segment <- function(
  x,
  metric = c("jaccard", "intersection", "ks", "mse", "chisq"),
  threshold = 0.5,
  stepsize = 1,
  max_sd_size = 1) {

  # Error checking
  if(!is.numeric(x)){
    stop("x must be a numeric vector")
  }
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
  # metric_func <- get(paste('histogram', metric, sep = "."))
  # TODO: pass in metric_func rather than metric as a param for fit_uniform?
  res <- lapply(seq(from = 1, to = num_bins - min_seg_size, by = stepsize), function(a) {
    lapply(seq(from = min_seg_size + a, to = num_bins, by = stepsize), function(b) {
      x_sub <- x[a:b]
      return(
        c(fit_uniform(x_sub, metric), list('seg_start' = a, 'seg_end' = b))
      )
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


# identify_mle_segment <- function(x) {
#   UseMethod('identify_mle_segment')
# }
#

# identify_mle_segment <- function(x) {
#   list(
#     "par" = NULL,
#     "dist" = "unif",
#     "metric" = "mle",
#     "value" = uniform.mle(metric, m),
#     "dens" = function(x = NULL, mpar = NULL, scale = TRUE) {
#       if(missing(x)) {
#         x <- bin
#       }
#       res <- ifelse(x >= min(bin) & x <= max(bin), p_unif[1], 0)
#       if(scale) res * N
#       else res
#     }
#   )
# }
