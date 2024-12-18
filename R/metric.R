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