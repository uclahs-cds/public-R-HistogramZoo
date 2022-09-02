
#' @export dgamma_flip
dgamma_flip <- function(x, shape, rate = 1, ...) {
  stats::dgamma(x = -x, shape = shape, rate = rate, ...)
}

#' @export pgamma_flip
pgamma_flip <- function(q, shape, rate = 1, lower.tail = TRUE, ...) {
  stats::pgamma(q = -q, shape = shape, rate = rate, lower.tail = !lower.tail, ...)
}

#' @export rgamma_flip
rgamma_flip <- function(n, shape, rate = 1, ...) {
  -stats::rgamma(n = n, shape = shape, rate = rate, ...)
}

#' @export qgamma_flip
qgamma_flip <- function(p, shape, rate = 1, lower.tail = TRUE, ...) {
  -stats::qgamma(p = p, shape = shape, rate = rate, lower.tail = !lower.tail, ...)
}
