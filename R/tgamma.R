dtgamma <- function(x, shape, rate = 1, a = 0, b = Inf, ...) {
  truncdist::dtrunc(x, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, ...)
}

ptgamma <- function(q, shape, rate = 1, a = 0, b = Inf, ...) {
  truncdist::ptrunc(q, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, ...)
}

rtgamma <- function(n, shape, rate = 1, a = 0, b = Inf, ...) {
  truncdist::rtrunc(n, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, ...)
}

qtgamma <- function(p, shape, rate = 1, a = 0, b = Inf, ...) {
  truncdist::qtrunc(p, spec = 'gamma', a = a, b = b, shape = shape, rate = rate, ...)
}
