dtgamma_flip <- function(x, shape, rate = 1, b, ...) {
  if(is.null(b)) stop("Need to supply upper bound.")
  truncdist::dtrunc(b - x, spec = 'gamma', a = 0, b = b, shape = shape, rate = rate, ...)
}

ptgamma_flip <- function(q, shape, rate = 1, b, lower.tail = TRUE, ...) {
  if(is.null(b)) stop("Need to supply upper bound.")
  truncdist::ptrunc(q, spec = 'gamma', a = 0, b = b, shape = shape, rate = rate, lower.tail = !lower.tail, ...)
}

rtgamma_flip <- function(n, shape, rate = 1, b, ...) {
  if(is.null(b)) stop("Need to supply upper bound.")
  b - truncdist::rtrunc(n, spec = 'gamma', a = 0, b = b, shape = shape, rate = rate, ...)
}

qtgamma_flip <- function(p, shape, rate = 1, b, ...) {
  if(is.null(b)) stop("Need to supply upper bound.")
  truncdist::qtrunc(p, spec = 'gamma', a = 0, b = b, shape = shape, rate = rate, ...)
}
