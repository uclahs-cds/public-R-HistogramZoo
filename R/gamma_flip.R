#' gamma flip
#'
#' @param x numeric vector of quantiles
#' @param shape shape parameter of gamma distribution
#' @param rate rate parameter of gamma distribution, default 1
#' @param offset offset of the right end of the distribution
#' @inheritParams stats::dgamma
#' @param ... additional parameters to be passed to stats::dgamma
#'
#' @return density of distribution
#' @export
dgamma_flip <- function(x, shape, rate = 1, offset = 0, ...) {
  stats::dgamma(x = offset - x, shape = shape, rate = rate, ...)
}

#' gamma flip
#'
#' @param q numeric vector of quantiles
#' @param shape shape parameter of gamma distribution
#' @param rate rate parameter of gamma distribution, default 1
#' @param offset offset of the right end of the distribution
#' @inheritParams stats::pgamma
#' @param ... additional parameters to be passed to stats::pgamma
#'
#' @return distribution function
#' @export
pgamma_flip <- function(q, shape, rate = 1, offset = 0, lower.tail = TRUE, ...) {
  stats::pgamma(q = offset - q, shape = shape, rate = rate, lower.tail = !lower.tail, ...)
}

#' gamma flip
#'
#' @param n number of observations
#' @param shape shape parameter of gamma distribution
#' @param rate rate parameter of gamma distribution, default 1
#' @param offset offset of the right end of the distribution
#' @inheritParams stats::rgamma
#' @param ... additional parameters to be passed to stats::rgamma
#'
#' @return random deviates
#' @export
rgamma_flip <- function(n, shape, rate = 1, offset = 0, ...) {
  offset - stats::rgamma(n = n, shape = shape, rate = rate, ...)
}

#' gamma flip
#'
#' @param p vector of probabilities
#' @param shape shape parameter of gamma distribution
#' @param rate rate parameter of gamma distribution, default 1
#' @param offset offset of the right end of the distribution
#' @inheritParams stats::qgamma
#' @param ... additional parameters to be passed to stats::qgamma
#'
#' @return quantile
#' @export
qgamma_flip <- function(p, shape, rate = 1, offset = 0, lower.tail = TRUE, ...) {
  offset - stats::qgamma(p = p, shape = shape, rate = rate, lower.tail = !lower.tail, ...)
}
