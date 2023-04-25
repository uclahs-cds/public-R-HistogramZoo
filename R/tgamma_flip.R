
#' The Truncated Flipped Gamma Distribution
#'
#' Density, distribution function, quantile function and random generation for the Truncated Flipped Gamma distribution,
#' a Gamma distribution reflected across the Y-axis, and truncted at endpoints [a, b]. Support provided by the `truncdist` package
#' with all parameters available for a standard Gamma distribution. One additional parameter, `offset` is provided to 
#' generate Truncated Flipped Gamma distributions with non-zero right endpoints, ending at the `offset` value.
#' 
#' @param x numeric vector of quantiles
#' @param q numeric vector of quantiles
#' @param n number of observations
#' @param p vector of probabilities
#' @param shape shape parameter of gamma distribution
#' @param rate rate parameter of gamma distribution, default 1
#' @param offset offset of the right end of the distribution
#' @param a distribution start
#' @param b distribution end
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param ... additional parameters to be passed to truncdist::dtrunc, truncdist::ptrunc, truncdist::rtrunc and truncdist::qtrunc
#' 
#' @return density of distribution
#' 
#' @rdname tgamma_flip
#' @export
dtgamma_flip <- function(x, shape, rate = 1, offset = 0, a = 0, b = Inf, ...) {
  truncdist::dtrunc(offset - x, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, ...)
}

#' @return distribution function
#' 
#' @rdname tgamma_flip
#' @export
ptgamma_flip <- function(q, shape, rate = 1, offset = 0, a = 0, b = Inf, lower.tail = TRUE, ...) {
  truncdist::ptrunc(offset - q, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, lower.tail = !lower.tail, ...)
}

#' @return random deviates
#' 
#' @rdname tgamma_flip
#' @export
rtgamma_flip <- function(n, shape, rate = 1, offset = 0, a = 0, b = Inf, ...) {
  offset - truncdist::rtrunc(n, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, ...)
}

#' @return quantile function
#' 
#' @rdname tgamma_flip
#' @export
qtgamma_flip <- function(p, shape, rate = 1, offset = 0, a = 0, b = Inf, lower.tail = TRUE, ...) {
  offset - truncdist::qtrunc(p, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, lower.tail = !lower.tail, ...)
}
