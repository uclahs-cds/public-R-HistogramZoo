#' Finds the largest uniform segment that is longer than threshold
#'
#' @param x numeric vector representing the density of a histogram
#' @param metric one of `jaccard`, `intersection`, `ks`, `mse`, `chisq` indicating metrics to use for fit optimization
#' @param threshold numeric, indicating the minimum proportion of the subsegment which should be tested
#' @param stepsize integer, indicating the stepsize (relative to the histogram bins) to take in the search for the uniform subsegment
#' @param max_sd_size numeric, the number of standard deviations of the computed metric distribution away from the optimal uniform which
#' has maximum length
#'
#' @return A data.frame with the following columns
#' \describe{
#'     \item{a}{start index of the maximum uniform segment}
#'     \item{b}{end index of the maximum uniform segment}
#'     \item{metric}{value of the fitted metric on the segment}
#'     \item{length}{length of the segment}
#'}
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
  metric_func <- get(paste('histogram', metric, sep = "."))

  p_unif <- generate_uniform_distribution(x)
  res <- lapply(seq(from = 1, to = num_bins - min_seg_size, by = stepsize), function(a) {
    lapply(seq(from = min_seg_size + a, to = num_bins, by = stepsize), function(b) {
      x_sub <- x[a:b]
      p_unif_sub <- generate_uniform_distribution(x_sub)
      h_sub <- x_sub / sum(x_sub)

      m <- metric_func(h_sub, p_unif_sub)

      return(
        list('start' = a, 'end' = b, 'metric' = m)
      )
    })
  })

  res_df <- do.call(rbind.data.frame, unlist(res, recursive = F))
  res_df$length <- res_df$end - res_df$start + 1

  # Select the longest interval that is within 1 sd of the maximum
  min_metric <- min(res_df$metric, na.rm = T)
  sd_metric <- stats::sd(res_df$metric, na.rm = T)
  sd_metric <- ifelse(is.na(sd_metric), 0, sd_metric) # If there's only 1 case
  # The range in which we are looking for the minimum
  res_sd_range <- res_df[res_df$metric <= min_metric + sd_metric * max_sd_size, ]
  max_interval_index <- which.max(res_sd_range$length)

  return(
    res_sd_range[max_interval_index,]
  )

}
