#' The Gamma Distribution
#' 
#' Density, distribution function, quantile function and random generation for the Gamma distribution
#' with `shape`, `rate` and other parameters inherited from the `stats` package. One additional parameter, `shift` 
#' is provided to generate Gamma distributions with non-zero starts, starting at the `shift` value.
#'
#' @param x numeric vector of quantiles
#' @param q numeric vector of quantiles
#' @param n number of observations
#' @param p vector of probabilities
#' @param shape shape parameter of gamma distribution
#' @param rate rate parameter of gamma distribution, default 1
#' @param shift shift of the left end of the distribution
#' @param ... additional parameters to be passed to stats::dgamma, stats::pgamma, stats::qgamma and stats::rgamma
#'
#' @inheritParams stats::dgamma
#' @return density of distribution
#' 
#' @rdname gamma
#' @export
dgamma <- function(x, shape, rate = 1, shift = 0, ...) {
  stats::dgamma(x = x - shift, shape = shape, rate = rate, ...)
}

#' @inheritParams stats::pgamma
#' @return distribution function
#' 
#' @rdname gamma
#' @export
pgamma <- function(q, shape, rate = 1, shift = 0, lower.tail = TRUE, ...) {
  stats::pgamma(q = q - shift, shape = shape, rate = rate, lower.tail = lower.tail, ...)
}

#' @inheritParams stats::rgamma
#' @return random deviates
#' 
#' @rdname gamma
#' @export
rgamma <- function(n, shape, rate = 1, shift = 0, ...) {
  shift + stats::rgamma(n = n, shape = shape, rate = rate, ...)
}

#' @inheritParams stats::qgamma
#' @return quantile
#' 
#' @rdname gamma
#' @export
qgamma <- function(p, shape, rate = 1, shift = 0, lower.tail = TRUE, ...) {
  shift + stats::qgamma(p = p, shape = shape, rate = rate, lower.tail = lower.tail, ...)
}
