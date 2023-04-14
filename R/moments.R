#' Computes the weighted sample standard deviation
#'
#' @export
#' @rdname weighted.sd
weighted.sd <- function(x, ...) UseMethod('weighted.sd')

#' Computes the weighted sample variance
#' @export
#' @rdname weighted.sd
weighted.var <- function(x, ...) UseMethod('weighted.var')

#' Computes the weighted sample skewness
#' @export
#' @rdname weighted.skewness
weighted.skewness <- function(x, ...) UseMethod('weighted.skewness')

#' Computes the sample moments of a Histogram or GenomicHistogram
#'
#' weighted.mean.Histogram computes a the sample variance of a Histogram or GenomicHistogram
#'
#' @param x Histogram or GenomicHistogram
#'
#' @exportS3Method weighted.mean Histogram
weighted.mean.Histogram <- function(x) {
  midpoints <- (x$interval_end + x$interval_start) / 2
  return(weighted.mean(x = midpoints, w = x$histogram_data))
}

#' weighted.var.Histogram computes a the sample variance of a Histogram or GenomicHistogram, with option for Sheppard's correction.
#'
#' @param x Histogram or GenomicHistogram
#' @param sheppard Should Sheppard's correction be applied? (subtract h^2/12 from the variance). Should only be applied when bin widths are all the same
#'
#' @rdname weighted.moments
#' @exportS3Method weighted.var Histogram
weighted.var.Histogram <- function(x, sheppard = FALSE) {
  midpoints <- (x$interval_end + x$interval_start) / 2
  adj <- 0
  if (sheppard) {
    if (length(unique(x$bin_width)) > 1 && sheppard) {
      stop('Non-equal bin widths with Sheppard\'s correction.')
    }
    adj <- (- x$bin_width[1]^2 / 12)
  }

  return(weighted.var(x = midpoints, w = x$histogram_data) + adj)
}

#' weighted.sd.Histogram computes a the sample variance of a Histogram or GenomicHistogram, with option for Sheppard's correction.
#'
#' @param x Histogram or GenomicHistogram
#'
#' @rdname weighted.moments
#' @exportS3Method weighted.sd Histogram
weighted.sd.Histogram <- function(x, sheppard = FALSE) {
  return(sqrt(weighted.var.Histogram(x, sheppard)))
}

#' Computes a the sample mean of a Histogram or GenomicHistogram
#'
#' @param x Histogram or GenomicHistogram
#'
#' @rdname weighted.moments
#' @exportS3Method weighted.skewness Histogram
weighted.skewness.Histogram <- function(x, ...) {
  midpoints <- (x$interval_end + x$interval_start) / 2
  weighted.skewness(x = midpoints, w = x$histogram_data, ...)
}

#' Computed the weighted sample standard deviation
#'
#' @param x values whose weighted sd is to be computed
#' @param w numerical vector of weight the same length as x
#' @export
weighted.sd.default <- function(x, w) {
  stopifnot(length(x) == length(w))
  return(sqrt(weighted.var(x, w)))
}

#' Computes the weighted sample skewness
#'
#' @param x values whose weighted variance is to be computed
#' @param w numerical vector of weight the same length as x
#' @export
weighted.var.default <- function(x, w, biased = FALSE) {
  stopifnot(length(x) == length(w))
  mu <- weighted.mean(x, w)
  N <- sum(w)
  denom <- if (biased) N else N - 1
  return((sum(w * x^2) - N * mu^2) / denom)
}

#' Computes the weighted sample skewness
#'
#' Either compute 'b' or 'g' skewness using the notation from: Joanes, D. N., & Gill, C. A. (1998). Comparing Measures of Sample Skewness and Kurtosis. Journal of the Royal Statistical Society.
#' 'b' skewness follows MINITAB and uses the unbiased formula for sd with (N - 1) in the denominator
#' 'g' skewness follows SAS format (and that of the R moments::skewness) and uses the biased formula for sd with N in the denominator
#' "For small samples from a normal distribution b has smaller mean-squared error than g... for small samples from non-normal distributions... g has a smaller mean-squared error."
#'
#' @param x values whose weighted skewness is to be computed
#' @param w numerical vector of weight the same length as x
#' @param type 'b' for (N - 1) adjustment to variance and 'g' for method of moments estimator
#'
#' @export
weighted.skewness.default <- function(x, w, type = c('b', 'g')) {
  stopifnot(length(x) == length(w))
  type <- match.arg(type)
  mu <- weighted.mean(x, w)
  N <- sum(w)
  third_central_moment <- sum(w * (x - mu)^3) / N
  biased <- type == 'g'
  return(third_central_moment / weighted.var(x, w, biased = biased)^1.5)
}
