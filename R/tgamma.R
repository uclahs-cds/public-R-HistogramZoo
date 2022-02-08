
#' truncated gamma
#'
#' @param x TODO
#' @param shape TODO
#' @param rate TODO
#' @param a TODO
#' @param b TODO
#' @param ... TODO
#'
#' @return
#' @export
dtgamma <- function(x, shape, rate = 1, a = 0, b = Inf, ...) {
  truncdist::dtrunc(x, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, ...)
}

#' truncated gamma
#'
#' @param q TODO
#' @param shape TODO
#' @param rate TODO
#' @param a TODO
#' @param b TODO
#' @param ... TODO
#'
#' @return
#' @export
ptgamma <- function(q, shape, rate = 1, a = 0, b = Inf, ...) {
  truncdist::ptrunc(q, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, ...)
}


#' truncated gamma
#'
#' @param n TODO
#' @param shape TODO
#' @param rate TODO
#' @param a TODO
#' @param b TODO
#' @param ... TODO
#'
#' @return
#' @export
rtgamma <- function(n, shape, rate = 1, a = 0, b = Inf, ...) {
  truncdist::rtrunc(n, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, ...)
}


#' truncated gamma
#'
#' @param p TODO
#' @param shape TODO
#' @param rate TODO
#' @param a TODO
#' @param b TODO
#' @param ... TODO
#'
#' @return
#' @export
qtgamma <- function(p, shape, rate = 1, a = 0, b = Inf, ...) {
  truncdist::qtrunc(p, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, ...)
}
