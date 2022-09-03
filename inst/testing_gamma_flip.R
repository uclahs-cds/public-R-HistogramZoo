

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

x = -c(1:20)
dgamma_flip_data = dgamma_flip(x, shape = 10)
plot(x, dgamma_flip_data)

dgamma_flip_data = dgamma_flip(rev(0:19), shape = 10, offset =20)
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
dtgamma_data = dtgamma(5:20, shape = 10, a = 5, b = 20)
plot(dtgamma_data)

dtgamma_flip_data = dtgamma_flip(5:20, shape = 10, a = 5, b = 20)
plot(dtgamma_flip_data)

# ptgamma
ptgamma_data = ptgamma(5:20, shape = 10, a = 5, b = 20)
plot(ptgamma_data)

x = rev(5:20)
ptgamma_flip_data = ptgamma_flip(rev(5:20), shape = 10, lower.tail = F, a = 5, b = 20)
plot(x, ptgamma_flip_data)

# qgamma
ptgamma(qtgamma(seq(0, 0.9, 0.1), shape = 10), shape = 10)
ptgamma(qtgamma(seq(0, 0.9, 0.1), shape = 10, a = 5, b = 20), shape = 10, a = 5, b = 20)
ptgamma_flip(qtgamma_flip(seq(0, 0.9, 0.1), shape = 10, a = 5, b = 20), shape = 10, a = 5, b = 20)

# Testing fit_distributions on gamma_flip ---------------------------------

set.seed(314)
histogram_data <- rgamma_flip(10000, shape = 10)
histogram <- observations_to_histogram(histogram_data)
histogram_data <- histogram$histogram_data

metric <- c("jaccard", "intersection", "ks", "mse", "chisq")
res <- fit_distributions(
  histogram_data,
  metric = metric,
  truncated = F,
  distributions = c("gamma_flip")
)
res_summary <- find_consensus_model(res)[['consensus']]

distribution_plotting_data <- lapply(res, function(m) {
  x <- seq(1, 26, by = 1)
  dens <- m$dens(x = seq_along(x), mpar = m$par)
  return(
    data.frame("fitted" = dens, "labels_x" = 1:26)
  )
})

res = segment_and_fit(histogram, distributions = c("gamma_flip"))
create_coverageplot(res)

histogram_data <- rtgamma_flip(10000, shape = 10, b = 20, a = 5)
histogram <- observations_to_histogram(histogram_data)
histogram_data <- histogram$histogram_data

metric <- c("jaccard", "intersection", "ks", "mse", "chisq")
res <- fit_distributions(
  histogram_data,
  metric = metric,
  truncated = T,
  distributions = c("gamma_flip")
)
res_summary <- find_consensus_model(res)[['consensus']]
distribution_plotting_data <- lapply(res, function(m) {
  x <- seq(1, 16, by = 1)
  dens <- m$dens(x = seq_along(x), mpar = m$par)
  return(
    data.frame("fitted" = dens, "labels_x" = 1:16)
  )
})

res = segment_and_fit(histogram, distributions = c("gamma_flip"), truncated = T, eps = 1, remove_low_entropy = F)
create_coverageplot(res)
