#' Histogram intersection
#' @param x histogram (vector of counts) 1
#' @param y histogram (vector of counts) 2
histogram.intersection <- function(x, y) {
  overlap = pmin(a = x, b = y, na.rm = T)
  # sums = c(sum(x), sum(y))
  # If one of the histograms is all zero, then return the other count
  1 - sum(overlap) / sum(x)
}

#' Jaccard index
#' @param x histogram (vector of counts) 1
#' @param y histogram (vector of counts) 2
histogram.jaccard <- function(x, y) {
  overlap = pmin(a = x, b = y, na.rm = T)
  union = pmax(a = x, b = y, na.rm = T)
  1 - sum(overlap)/sum(union)
}

#' Kolmogorov-Smirnov Divergance
#' @param x histogram (vector of counts) 1
#' @param y histogram (vector of counts) 2
histogram.ks <- function(x, y) {
  max(abs(x - y), na.rm = TRUE)
}

#' Mean-squared error
#' @param x histogram (vector of counts) 1
#' @param y histogram (vector of counts) 2
histogram.mse <- function(x, y) {
  mean((x - y)^2)
}

#' Chi-Squared
#' @param x histogram (vector of counts) 1
#' @param y histogram (vector of counts) 2
histogram.chisq <- function(x, y) {
  sum((x - y)^2 / (x + y))
}
