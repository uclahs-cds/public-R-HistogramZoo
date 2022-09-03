
#' truncated gamma flip
#'
#' @param x numeric vector of quantiles
#' @param shape shape parameter of gamma distribution
#' @param rate rate parameter of gamma distribution, default 1
#' @param a distribution start
#' @param b distribution end
#' @param ... additional parameters to be passed to truncdist::dtrunc
#'
#' @return density of distribution
#' @export
dtgamma_flip <- function(x, shape, rate = 1, a = 0, b, ...) {
  if(is.null(b)) stop("Need to supply upper bound.")
  truncdist::dtrunc(a + b - x, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, ...)
}

#' truncated gamma flip
#'
#' @param q numeric vector of quantiles
#' @param shape shape parameter of gamma distribution
#' @param rate rate parameter of gamma distribution, default 1
#' @param a distribution start
#' @param b distribution end
#' @param ... additional parameters to be passed to truncdist::dtrunc
#'
#' @return distribution function
#' @export
ptgamma_flip <- function(q, shape, rate = 1, a = 0, b, lower.tail = TRUE, ...) {
  if(is.null(b)) stop("Need to supply upper bound.")
  truncdist::ptrunc(a + b - q, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, lower.tail = !lower.tail, ...)
}

#' truncated gamma flip
#'
#' @param n number of observations
#' @param shape shape parameter of gamma distribution
#' @param rate rate parameter of gamma distribution, default 1
#' @param a distribution start
#' @param b distribution end
#' @param ... additional parameters to be passed to truncdist::dtrunc
#'
#' @return random deviates
#' @export
rtgamma_flip <- function(n, shape, rate = 1, a = 0, b, ...) {
  if(is.null(b)) stop("Need to supply upper bound.")
  a + b - truncdist::rtrunc(n, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, ...)
}

#' truncated gamma
#'
#' @param p vector of probabilities
#' @param shape shape parameter of gamma distribution
#' @param rate rate parameter of gamma distribution, default 1
#' @param a distribution start
#' @param b distribution end
#' @param ... additional parameters to be passed to truncdist::dtrunc
#'
#' @return quantile function
#' @export
qtgamma_flip <- function(p, shape, rate = 1, a = 0, b, lower.tail = TRUE, ...) {
  if(is.null(b)) stop("Need to supply upper bound.")
  a + b - truncdist::qtrunc(p, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, lower.tail = !lower.tail, ...)
}
