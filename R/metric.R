#' Histogram intersection
#'
#' @param x histogram (vector of counts) 1
#' @param y histogram (vector of counts) 2
#'
#' @return the intersection of x and y
histogram.intersection <- function(x, y) {
  overlap = pmin(a = x, b = y, na.rm = T)
  # sums = c(sum(x), sum(y))
  # If one of the histograms is all zero, then return the other count
  1 - sum(overlap) / sum(x)
}

#' Jaccard index
#'
#' @param x histogram (vector of counts) 1
#' @param y histogram (vector of counts) 2
#'
#' @return the jaccard similarity of x and y
histogram.jaccard <- function(x, y) {
  overlap = pmin(a = x, b = y, na.rm = T)
  union = pmax(a = x, b = y, na.rm = T)
  1 - sum(overlap)/sum(union)
}

#' Kolmogorov-Smirnov divergence
#' @param x histogram (vector of counts) 1
#' @param y histogram (vector of counts) 2
#'
#' @return the Kolmogorov-Smirnov divergence between x and y
histogram.ks <- function(x, y) {
  max(abs(x - y), na.rm = TRUE)
}

#' Mean-squared error
#'
#' @param x histogram (vector of counts) 1
#' @param y histogram (vector of counts) 2
#'
#' @return the mean-squared error of x and y
histogram.mse <- function(x, y) {
  mean((x - y)^2)
}

#' Chi-squared
#'
#' @param x histogram (vector of counts) 1
#' @param y histogram (vector of counts) 2
#'
#' @return the chi-squared distance between x and y
histogram.chisq <- function(x, y) {
  sum((x - y)^2 / (x + y))
}

uniform_mle <- function(x, a, b, inclusive = TRUE, log = TRUE) {
  UseMethod('uniform_mle')
}

# Internal, generalizes functionality between uniform_mle methods
uniform_mle_helper <- function(x, x.start, x.end, a, b, inclusive = TRUE, log = TRUE) {
  N <- sum(x)
  if (inclusive && any(x.start < a | x.end > b)) {
    return(-Inf)
  } else if (! inclusive && any(x.start <= a & x.end >= b)) {
    return(-Inf)
  }

  if (log) {
    return(- N * log(b - a))
  } else {
    return((1 / (b - a))^N)
  }
}

# x numeric
# Assumes bins are (0, 1], (1, 2], ... (N - 1, N]
#' @exportS3Method uniform_mle numeric
uniform_mle.numeric <- function(x, a, b, inclusive = TRUE, log = TRUE) {
  L <- length(x)
  x.end <- 1:L
  x.start <- x.end - 1

  return(uniform_mle_helper(x, x.start, x.end, a, b, inclusive, log))
}

# x Histogram
#' @exportS3Method uniform_mle Histogram
uniform_mle.Histogram <- function(x, a, b, inclusive = TRUE, log = TRUE) {
  x.start <- head(x$interval_start, n = 1)
  x.end <- tail(x$interval_end, n = 1)

  return(uniform_mle_helper(x$histogram_data, x.start, x.end, a, b, inclusive, log))
}
