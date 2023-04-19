
#' truncated gamma
#'
#' @param x numeric vector of quantiles
#' @param shape shape parameter of gamma distribution
#' @param rate rate parameter of gamma distribution, default 1
#' @param shift shift of the left end of the distribution
#' @param a distribution start
#' @param b distribution end
#' @param ... additional parameters to be passed to truncdist::dtrunc
#'
#' @return density of distribution
#' @export
dtgamma <- function(x, shape, rate = 1, shift = 0, a = 0, b = Inf, ...) {
  truncdist::dtrunc(x - shift, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, ...)
}

#' truncated gamma
#'
#' @param q numeric vector of quantiles
#' @param shape shape parameter of gamma distribution
#' @param rate rate parameter of gamma distribution, default 1
#' @param shift shift of the left end of the distribution
#' @param a distribution start
#' @param b distribution end
#' @param ... additional parameters to be passed to truncdist::dtrunc
#'
#' @return distribution function
#' @export
ptgamma <- function(q, shape, rate = 1, shift = 0, a = 0, b = Inf, ...) {
  truncdist::ptrunc(q - shift, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, ...)
}

#' truncated gamma
#'
#' @param n number of observations
#' @param shape shape parameter of gamma distribution
#' @param rate rate parameter of gamma distribution, default 1
#' @param shift shift of the left end of the distribution
#' @param a distribution start
#' @param b distribution end
#' @param ... additional parameters to be passed to truncdist::dtrunc
#'
#' @return random deviates
#' @export
rtgamma <- function(n, shape, rate = 1, shift = 0, a = 0, b = Inf, ...) {
  shift + truncdist::rtrunc(n, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, ...)
}

#' truncated gamma
#'
#' @param p vector of probabilities
#' @param shape shape parameter of gamma distribution
#' @param rate rate parameter of gamma distribution, default 1
#' @param shift shift of the left end of the distribution
#' @param a distribution start
#' @param b distribution end
#' @param ... additional parameters to be passed to truncdist::dtrunc
#'
#' @return quantile function
#' @export
qtgamma <- function(p, shape, rate = 1, shift = 0, a = 0, b = Inf, ...) {
  shift + truncdist::qtrunc(p, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, ...)
}
