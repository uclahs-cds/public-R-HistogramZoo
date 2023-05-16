#' The Flipped Gamma Distribution
#' 
#' Density, distribution function, quantile function and random generation for the Flipped Gamma distribution,
#' a Gamma distribution reflected across the Y-axis. `shape`, `rate` and other parameters are inherited from the `stats` package. 
#' One additional parameter, `offset` is provided to generate Flipped Gamma distributions with non-zero right endpoints, ending
#' at the `offset` value.
#'
#' @param x numeric vector of quantiles
#' @param q numeric vector of quantiles
#' @param n number of observations
#' @param p vector of probabilities
#' @param shape shape parameter of gamma distribution
#' @param rate rate parameter of gamma distribution, default 1
#' @param offset offset of the right end of the distribution
#' @param ... additional parameters to be passed to stats::dgamma, stats::pgamma, stats::rgamma and stats::qgamma
#'
#' @inheritParams stats::dgamma
#' @return density of distribution
#' 
#' @rdname gamma_flip
#' @export
dgamma_flip <- function(x, shape, rate = 1, offset = 0, ...) {
  stats::dgamma(x = offset - x, shape = shape, rate = rate, ...)
}

#' @inheritParams stats::pgamma
#' @return distribution function
#'
#' @rdname gamma_flip
#' @export
pgamma_flip <- function(q, shape, rate = 1, offset = 0, lower.tail = TRUE, ...) {
  stats::pgamma(q = offset - q, shape = shape, rate = rate, lower.tail = !lower.tail, ...)
}


#' @inheritParams stats::rgamma
#' @return random deviates
#'
#' @rdname gamma_flip
#' @export
rgamma_flip <- function(n, shape, rate = 1, offset = 0, ...) {
  offset - stats::rgamma(n = n, shape = shape, rate = rate, ...)
}

#' @inheritParams stats::qgamma
#' @return quantile
#'
#' @rdname gamma_flip
#' @export
qgamma_flip <- function(p, shape, rate = 1, offset = 0, lower.tail = TRUE, ...) {
  offset - stats::qgamma(p = p, shape = shape, rate = rate, lower.tail = !lower.tail, ...)
}
