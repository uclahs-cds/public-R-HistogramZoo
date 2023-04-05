context("fit_distributions_MLE")

# Initializing
metric <- c("jaccard", "intersection", "ks", "mse", "chisq")
distributions <- c("norm", "gamma", "unif", "gamma_flip")

test_that("fit_distributions: MLE - normal", {
  set.seed(314)
  histogram_data <- rnorm(10000, mean = 0, sd = 5)
  histogram <- observations_to_histogram(histogram_data)
  histogram_data <- histogram$histogram_data
  midpoint <- length(histogram_data)/2

  res <- fit_distributions(
    histogram_data,
    metric = 'mle',
    truncated = F,
    distributions = distributions
  )

  res_summary <- find_consensus_model(res, metric = 'mle')[['mle']]

  expect_equal(res_summary$dist, "norm")
  expect_true(res_summary$par$mean < (midpoint + 1) & res_summary$par$mean > (midpoint - 1))
  expect_true(res_summary$par$sd < (5 + 1) & res_summary$par$sd > (5 - 1))
  expect_true(res_summary$value > 0.8)
})

test_that("fit_distributions: MLE - gamma", {
  set.seed(314)
  shape <- 2
  rate <- 0.1
  histogram_data <- rgamma(10000, shape=shape, rate=rate)
  histogram <- observations_to_histogram(histogram_data)
  histogram_data <- histogram$histogram_data

  res <- fit_distributions(
    histogram_data,
    metric = 'mle',
    truncated = F,
    distributions = distributions
  )

  res_summary <- find_consensus_model(res, metric = 'mle')[['mle']]

  expect_equal(res_summary$dist, "gamma")
  expect_true(res_summary$par$rate < (rate + 0.05) & res_summary$par$rate > (rate - 0.05))
  expect_true(res_summary$par$shape < (shape + 1) & res_summary$par$shape > (shape - 1))
  expect_true(res_summary$value > 0.8)
})

test_that("fit_distributions: MLE - unif", {

  set.seed(314)
  histogram_data <- round(runif(10, min = 0, max = 1)*100) + rep(100, 10)

  res <- fit_distributions(
    histogram_data,
    metric = 'mle',
    truncated = F,
    distributions = distributions
  )

  res_summary <- find_consensus_model(res, metric = 'mle')[['mle']]

  expect_equal(res_summary$dist, "unif")
  expect_true(res_summary$value > 0.8)

})
