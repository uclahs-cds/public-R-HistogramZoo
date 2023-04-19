
# dgamma: density function of the gamma distribution
# pgamma: cumulative density function of the gamma distribution
# qgamma: quantile function of the gamma distribution
# rgamma: random sampling from the gamma distribution

dgamma <- function(x, shape, rate = 1, shift = 0, ...) {
  stats::dgamma(x = x - shift, shape = shape, rate = rate, ...)
}

pgamma <- function(q, shape, rate = 1, shift = 0, lower.tail = TRUE, ...) {
  stats::pgamma(q = q - shift, shape = shape, rate = rate, lower.tail = lower.tail, ...)
}

rgamma <- function(n, shape, rate = 1, shift = 0, ...) {
  shift + stats::rgamma(n = n, shape = shape, rate = rate, ...)
}

qgamma <- function(p, shape, rate = 1, shift = 0, lower.tail = TRUE, ...) {
  shift + stats::qgamma(p = p, shape = shape, rate = rate, lower.tail = lower.tail, ...)
}

x <- dgamma(seq(0, 30), shape = 10, shift = 10)
plot(x)

x <- pgamma(seq(1, 100), shape = 10, shift = 50)
plot(x)

x <- rgamma(10000, shape = 10, shift = 10)
hist(x)

x <- qgamma(seq(0, 1, 0.1), shape = 10, shift = 10)
plot(x)