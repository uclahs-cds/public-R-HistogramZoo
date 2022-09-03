

library(HistogramZoo)

# Preamble ----------------------------------------------------------------
# Testing gamma_flip
# Also testing truncated_gamma_flip and truncated_gamma

# Testing -----------------------------------------------------------------

set.seed(314)

# rgamma
gamma_data = rgamma(10000, shape = 10)
hist(gamma_data, breaks = 100)

rgamma_data = rgamma_flip(10000, shape = 10)
hist(rgamma_data, breaks = 100)

rgamma_data = rgamma_flip(10000, shape = 10, offset = 20)
hist(rgamma_data, breaks = 100)

# dgamma
dgamma_data = dgamma(1:20, shape = 10)
plot(dgamma_data)

x = -rev(1:20)
dgamma_flip_data = dgamma_flip(x, shape = 10)
plot(x, dgamma_flip_data)

dgamma_flip_data = dgamma_flip(1:20, shape = 10, offset =20)
plot(dgamma_flip_data)

# pgamma
pgamma_data = pgamma(1:20, shape = 10)
plot(pgamma_data)

x = -c(1:20)
pgamma_flip_data = pgamma_flip(x, shape = 10, lower.tail = F)
plot(x, pgamma_flip_data)

x = -c(1:20)+20
pgamma_flip_data = pgamma_flip(x, shape = 10, offset = 20, lower.tail = F)
plot(x, pgamma_flip_data)

# qgamma
pgamma(qgamma(seq(0, 0.9, 0.1), shape = 10), shape = 10)
pgamma_flip(qgamma_flip(seq(0, 0.9, 0.1), shape = 10), shape = 10)
pgamma_flip(qgamma_flip(seq(0, 0.9, 0.1), shape = 10, offset = 20), shape = 10, offset = 20)

# Testing truncated gamma -------------------------------------------------

# rtgamma
rtgamma_data = rtgamma(10000, shape = 10, a = 5, b = 20)
hist(rtgamma_data, breaks = 100)

rtgamma_flip_data = rtgamma_flip(10000, shape = 10, a = 5, b = 20)
hist(rtgamma_flip_data, breaks = 100)

# dtgamma
dtgamma_data = dtgamma(1:25, shape = 10, a = 5, b = 20)
plot(dtgamma_data)

dtgamma_flip_data = dtgamma_flip(1:25, shape = 10, a = 5, b = 20)
plot(dtgamma_flip_data)

# ptgamma
ptgamma_data = ptgamma(5:20, shape = 10, a = 5, b = 20)
plot(ptgamma_data)

x = -c(0:15)+20
ptgamma_flip_data = ptgamma_flip(x, shape = 10, lower.tail = F, a = 5, b = 20)
plot(x, ptgamma_flip_data)

# qgamma
ptgamma(qtgamma(seq(0, 0.9, 0.1), shape = 10), shape = 10)
ptgamma(qtgamma(seq(0, 0.9, 0.1), shape = 10, a = 5, b = 20), shape = 10, a = 5, b = 20)
ptgamma_flip(qtgamma_flip(seq(0, 0.9, 0.1), shape = 10, a = 5, b = 20), shape = 10, a = 5, b = 20)
