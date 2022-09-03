context("fit_distributions")


# Initializing
metric <- c("jaccard", "intersection", "ks", "mse", "chisq")
distributions <- c("norm", "gamma", "gamma_flip", "unif")
fit_names <- c("par", "dist", "metric", "value", "dens")

test_that("fit_distributions works ", {
  
  histogram_data <- round(rnorm(100, mean = 0, sd = 5))
  histogram_data <- table(histogram_data)

  res <- fit_distributions(
    histogram_data,
    metric = metric,
    truncated = F,
    distributions = distributions
  )

  expect_length(res, 20)
  expect_true(
    all(distributions %in% unlist(lapply(res, `[`, "dist")))
  )
  expect_true(
    all(metric %in% unlist(lapply(res, `[`, "metric")))
  )
  expect_named(res[[1]], fit_names)

})

test_that("fit_distributions: normal", {
  
  set.seed(314)
  histogram_data <- rnorm(10000, mean = 0, sd = 5)
  histogram <- observations_to_histogram(histogram_data)
  histogram_data <- histogram$histogram_data
  midpoint <- length(histogram_data)/2
  
  res <- fit_distributions(
    histogram_data,
    metric = metric,
    truncated = F,
    distributions = distributions
  )
  
  res_summary <- find_consensus_model(res)[['consensus']]
  
  expect_equal(res_summary$dist, "norm")
  expect_true(res_summary$par$mean < (midpoint + 1) & res_summary$par$mean > (midpoint - 1))
  expect_true(res_summary$par$sd < (5 + 1) & res_summary$par$sd > (5 - 1))
  expect_true(res_summary$value > 0.8)
  
  # Truncated
  histogram_data <- histogram_data[5:length(histogram_data)-5]
  
  res <- fit_distributions(
    histogram_data,
    metric = metric,
    truncated = T,
    distributions = distributions
  )
  
  expect_equal(res_summary$dist, "norm")
  expect_true(res_summary$par$mean < (midpoint + 1) & res_summary$par$mean > (midpoint - 1))
  expect_true(res_summary$par$sd < (5 + 1) & res_summary$par$sd > (5 - 1))
  expect_true(res_summary$value > 0.8)
  
})

test_that("fit_distributions: unif", {
  
  set.seed(314)
  histogram_data <- round(runif(10, min = 0, max = 1)*100) + rep(100, 10)
  
  res <- fit_distributions(
    histogram_data,
    metric = metric,
    truncated = F,
    distributions = distributions
  )
  
  res_summary <- find_consensus_model(res)[['consensus']]
  
  expect_equal(res_summary$dist, "unif")
  expect_true(res_summary$value > 0.8)
  
})

test_that("fit_distributions: gamma", {
  
  set.seed(314)
  shape <- 2
  rate <- 0.1
  histogram_data <- rgamma(10000, shape=shape, rate=rate)
  histogram <- observations_to_histogram(histogram_data)
  histogram_data <- histogram$histogram_data
  
  res <- fit_distributions(
    histogram_data,
    metric = metric,
    truncated = F,
    distributions = distributions
  )
  
  res_summary <- find_consensus_model(res)[['consensus']]
  
  expect_equal(res_summary$dist, "gamma")
  expect_true(res_summary$par$rate < (rate + 0.05) & res_summary$par$rate > (rate - 0.05))
  expect_true(res_summary$par$shape < (shape + 1) & res_summary$par$shape > (shape - 1))
  expect_true(res_summary$value > 0.8)
  
  # Truncated
  histogram_data <- histogram_data[10:length(histogram_data)-10]
  
  res <- fit_distributions(
    histogram_data,
    metric = metric,
    truncated = T,
    distributions = distributions
  )
  
  res_summary <- find_consensus_model(res)[['consensus']]
  
  expect_equal(res_summary$dist, "gamma")
  expect_true(res_summary$par$rate < (rate + 0.05) & res_summary$par$rate > (rate - 0.05))
  expect_true(res_summary$par$shape < (shape + 1) & res_summary$par$shape > (shape - 1))
  expect_true(res_summary$value > 0.8)
  
})

test_that("fit_distributions: gamma_flip", {
  
  set.seed(314)
  shape <- 2
  rate <- 0.1
  histogram_data <- rgamma_flip(10000, shape=shape, rate=rate, offset = 136)
  histogram <- observations_to_histogram(histogram_data)
  histogram_data <- histogram$histogram_data
  
  res <- fit_distributions(
    histogram_data,
    metric = metric,
    truncated = F,
    distributions = distributions
  )
  
  res_summary <- find_consensus_model(res)[['consensus']]
  
  expect_equal(res_summary$dist, "gamma_flip")
  expect_true(res_summary$par$rate < (rate + 0.05) & res_summary$par$rate > (rate - 0.05))
  expect_true(res_summary$par$shape < (shape + 1) & res_summary$par$shape > (shape - 1))
  expect_true(res_summary$value > 0.8)
  
  # Truncated
  histogram_data <- histogram_data[10:length(histogram_data)-10]
  
  res <- fit_distributions(
    histogram_data,
    metric = metric,
    truncated = T,
    distributions = distributions
  )
  
  res_summary <- find_consensus_model(res)[['consensus']]
  
  expect_equal(res_summary$dist, "gamma_flip")
  expect_true(res_summary$par$rate < (rate + 0.05) & res_summary$par$rate > (rate - 0.05))
  expect_true(res_summary$par$shape < (shape + 1) & res_summary$par$shape > (shape - 1))
  expect_true(res_summary$value > 0.8)
  
})

