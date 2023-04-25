#' gamma
#'
#' @param x numeric vector of quantiles
#' @param shape shape parameter of gamma distribution
#' @param rate rate parameter of gamma distribution, default 1
#' @param shift shift of the left end of the distribution
#' @param ... additional parameters to be passed to stats::dgamma
#'
#' @return density of distribution
#' @export
dgamma <- function(x, shape, rate = 1, shift = 0, ...) {
  stats::dgamma(x = x - shift, shape = shape, rate = rate, ...)
}

#' gamma
#'
#' @param q numeric vector of quantiles
#' @param shape shape parameter of gamma distribution
#' @param rate rate parameter of gamma distribution, default 1
#' @param shift shift of the left end of the distribution
#' @inheritParams stats::pgamma
#' @param ... additional parameters to be passed to stats::pgamma
#'
#' @return distribution function
#' @export
pgamma <- function(q, shape, rate = 1, shift = 0, lower.tail = TRUE, ...) {
  stats::pgamma(q = q - shift, shape = shape, rate = rate, lower.tail = lower.tail, ...)
}


#' gamma
#'
#' @param n number of observations
#' @param shape shape parameter of gamma distribution
#' @param rate rate parameter of gamma distribution, default 1
#' @param shift shift of the left end of the distribution
#' @param ... additional parameters to be passed to stats::rgamma
#'
#' @return random deviates
#' @export
rgamma <- function(n, shape, rate = 1, shift = 0, ...) {
  shift + stats::rgamma(n = n, shape = shape, rate = rate, ...)
}

#' gamma
#'
#' @param p vector of probabilities
#' @param shape shape parameter of gamma distribution
#' @param rate rate parameter of gamma distribution, default 1
#' @param shift shift of the left end of the distribution
#' @inheritParams stats::qgamma
#' @param ... additional parameters to be passed to stats::qgamma
#'
#' @return quantile
#' @export
qgamma <- function(p, shape, rate = 1, shift = 0, lower.tail = TRUE, ...) {
  shift + stats::qgamma(p = p, shape = shape, rate = rate, lower.tail = lower.tail, ...)
}
