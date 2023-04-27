
#' The Truncated Gamma Distribution
#' 
#' Density, distribution function, quantile function and random generation for the Truncated Gamma distribution with endpoints
#' at [a, b]. Support provided by the `truncdist` package with all parameters available for a standard Gamma distribution. 
#' One additional parameter, `shift` is provided to generate Gamma distributions with non-zero starts, starting at the `shift` value.
#'
#' @param x numeric vector of quantiles
#' @param q numeric vector of quantiles
#' @param n number of observations
#' @param p vector of probabilities
#' @param shape shape parameter of gamma distribution
#' @param rate rate parameter of gamma distribution, default 1
#' @param shift shift of the left end of the distribution
#' @param a distribution start
#' @param b distribution end
#' @param ... additional parameters to be passed to truncdist::dtrunc, truncdist::ptrunc, truncdist::rtrunc and truncdist::qtrunc
#'
#' @return density of distribution
#' 
#' @rdname tgamma
#' @export
dtgamma <- function(x, shape, rate = 1, shift = 0, a = 0, b = Inf, ...) {
  truncdist::dtrunc(x - shift, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, ...)
}

#' @return distribution function
#'
#' @rdname tgamma
#' @export
ptgamma <- function(q, shape, rate = 1, shift = 0, a = 0, b = Inf, ...) {
  truncdist::ptrunc(q - shift, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, ...)
}

#' @return random deviates
#'
#' @rdname tgamma
#' @export
rtgamma <- function(n, shape, rate = 1, shift = 0, a = 0, b = Inf, ...) {
  shift + truncdist::rtrunc(n, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, ...)
}

#' @return quantile function
#'
#' @rdname tgamma
#' @export
qtgamma <- function(p, shape, rate = 1, shift = 0, a = 0, b = Inf, ...) {
  shift + truncdist::qtrunc(p, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, ...)
}
