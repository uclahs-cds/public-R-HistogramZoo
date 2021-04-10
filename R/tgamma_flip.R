
#' @export dtgamma_flip
dtgamma_flip <- function(x, shape, rate = 1, b, ...) {
  if(is.null(b)) stop("Need to supply upper bound.")
  truncdist::dtrunc(b - x, spec = 'gamma', a = 0, b = b, shape = shape, rate = rate, ...)
}

#' @export ptgamma_flip
ptgamma_flip <- function(q, shape, rate = 1, b, lower.tail = TRUE, ...) {
  if(is.null(b)) stop("Need to supply upper bound.")
  truncdist::ptrunc(q, spec = 'gamma', a = 0, b = b, shape = shape, rate = rate, lower.tail = !lower.tail, ...)
}

#' @export rtgamma_flip
rtgamma_flip <- function(n, shape, rate = 1, b, ...) {
  if(is.null(b)) stop("Need to supply upper bound.")
  b - truncdist::rtrunc(n, spec = 'gamma', a = 0, b = b, shape = shape, rate = rate, ...)
}

#' @export qtgamma_flip
qtgamma_flip <- function(p, shape, rate = 1, b, ...) {
  if(is.null(b)) stop("Need to supply upper bound.")
  truncdist::qtrunc(p, spec = 'gamma', a = 0, b = b, shape = shape, rate = rate, ...)
}
