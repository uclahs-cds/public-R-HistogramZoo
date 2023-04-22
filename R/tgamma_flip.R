
#' truncated gamma flip
#'
#' @param x numeric vector of quantiles
#' @param shape shape parameter of gamma distribution
#' @param rate rate parameter of gamma distribution, default 1
#' @param offset offset of the right end of the distribution
#' @param a distribution start
#' @param b distribution end
#' @param ... additional parameters to be passed to truncdist::dtrunc
#'
#' @return density of distribution
#' @export
dtgamma_flip <- function(x, shape, rate = 1, offset = 0, a = 0, b = Inf, ...) {
  truncdist::dtrunc(offset - x, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, ...)
}

#' truncated gamma flip
#'
#' @param q numeric vector of quantiles
#' @param shape shape parameter of gamma distribution
#' @param rate rate parameter of gamma distribution, default 1
#' @param offset offset of the right end of the distribution
#' @param a distribution start
#' @param b distribution end
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param ... additional parameters to be passed to truncdist::dtrunc
#'
#' @return distribution function
#' @export
ptgamma_flip <- function(q, shape, rate = 1, offset = 0, a = 0, b = Inf, lower.tail = TRUE, ...) {
  truncdist::ptrunc(offset - q, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, lower.tail = !lower.tail, ...)
}

#' truncated gamma flip
#'
#' @param n number of observations
#' @param shape shape parameter of gamma distribution
#' @param rate rate parameter of gamma distribution, default 1
#' @param offset offset of the right end of the distribution
#' @param a distribution start
#' @param b distribution end
#' @param ... additional parameters to be passed to truncdist::dtrunc
#'
#' @return random deviates
#' @export
rtgamma_flip <- function(n, shape, rate = 1, offset = 0, a = 0, b = Inf, ...) {
  offset - truncdist::rtrunc(n, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, ...)
}

#' truncated gamma
#'
#' @param p vector of probabilities
#' @param shape shape parameter of gamma distribution
#' @param rate rate parameter of gamma distribution, default 1
#' @param offset offset of the right end of the distribution
#' @param a distribution start
#' @param b distribution end
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
#' @param ... additional parameters to be passed to truncdist::dtrunc
#'
#' @return quantile function
#' @export
qtgamma_flip <- function(p, shape, rate = 1, offset = 0, a = 0, b = Inf, lower.tail = TRUE, ...) {
  offset - truncdist::qtrunc(p, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, lower.tail = !lower.tail, ...)
}
