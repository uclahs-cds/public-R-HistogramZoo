#' Weighted higher order moments
#'
#' Computes weighted sample mean, standard deviation, variance, and skewness.
#' See stats::\link[stats]{weighted.mean}.
#'
#' @export
#' @rdname weighted.moments
weighted.sd <- function(x, ...) UseMethod('weighted.sd')

#' @export
#' @rdname weighted.moments
weighted.var <- function(x, ...) UseMethod('weighted.var')

#' @export
#' @rdname weighted.moments
weighted.skewness <- function(x, ...) UseMethod('weighted.skewness')

#' @param x Histogram or GenomicHistogram
#'
#' @exportS3Method weighted.mean Histogram
#' @rdname weighted.moments
weighted.mean.Histogram <- function(x) {
  midpoints <- find_midpoint(x$interval_start, x$interval_end)
  return(weighted.mean(x = midpoints, w = x$histogram_data))
}

#' Sample variance of a Histogram or GenomicHistogram
#'
#' @param x Histogram or GenomicHistogram
#' @param sheppard Should Sheppard's correction be applied? (subtract h^2/12 from the variance). Should only be applied when bin widths are all the same
#'
#' @rdname weighted.moments
#' @exportS3Method weighted.var Histogram
weighted.var.Histogram <- function(x, sheppard = FALSE) {
  bin_width <- x$interval_end - x$interval_start
  midpoints <- find_midpoint(x$interval_start, x$interval_end)
  adj <- 0
  if (sheppard) {
    if (length(unique(bin_width)) > 1 && sheppard) {
      stop('Non-equal bin widths with Sheppard\'s correction.')
    }
    adj <- (- bin_width[1]^2 / 12)
  }

  return(weighted.var(x = midpoints, w = x$histogram_data) + adj)
}

#' Sample sd of a Histogram or GenomicHistogram
#'
#' @param x Histogram or GenomicHistogram
#'
#'
#' @rdname weighted.moments
#' @exportS3Method weighted.sd Histogram
weighted.sd.Histogram <- function(x, sheppard = FALSE) {
  return(sqrt(weighted.var.Histogram(x, sheppard)))
}

#' Sample skewness of a Histogram or GenomicHistogram
#'
#' @param x Histogram or GenomicHistogram
#'
#' @rdname weighted.moments
#' @exportS3Method weighted.skewness Histogram
#' @examples
#' x <- Histogram(c(1,1,2,2,3,3,4,3,2,1))
#' weighted.mean(x)
#' weighted.sd(x)
#' weighted.var(x)
#' weighted.skewness(x)
weighted.skewness.Histogram <- function(x, ...) {
  midpoints <- find_midpoint(x$interval_start, x$interval_end)
  weighted.skewness(x = midpoints, w = x$histogram_data, ...)
}

#' @param x values whose weighted sd is to be computed
#' @param w numerical vector of weight the same length as x
#' @export
#' @rdname weighted.moments
weighted.sd.default <- function(x, w) {
  stopifnot(length(x) == length(w))
  return(sqrt(weighted.var(x, w)))
}

#' @param x values whose weighted variance is to be computed
#' @param w numerical vector of weight the same length as x
#' @export
#' @rdname weighted.moments
weighted.var.default <- function(x, w, biased = FALSE) {
  stopifnot(length(x) == length(w))
  mu <- weighted.mean(x, w)
  N <- sum(w)
  denom <- if (biased) N else N - 1
  return((sum(w * x^2) - N * mu^2) / denom)
}

#' @section Sample Skewness:
#' Sample skewness either computes 'b' or 'g' skewness using the notation from [1].
#' \itemize{
#'   \item `b` skewness follows MINITAB and uses the unbiased formula for sd with (N - 1) in the denominator
#'   \item `g` skewness follows SAS format (and that of the R moments::\link[moments]{skewness}) and uses the biased formula for sd with N in the denominator
#' }
#' "For small samples from a normal distribution b has smaller mean-squared error than g... for small samples from non-normal distributions... g has a smaller mean-squared error."
#'
#' @references
#' [1] Joanes, D. N., & Gill, C. A. (1998). Comparing Measures of Sample Skewness and Kurtosis. Journal of the Royal Statistical Society. Series D (The Statistician), 47(1), 183–189. http://www.jstor.org/stable/2988433
#'
#' @param x values whose weighted skewness is to be computed
#' @param w numerical vector of weight the same length as x
#' @param type 'b' for (N - 1) adjustment to variance and 'g' for method of moments estimator.
#' @rdname weighted.moments
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
