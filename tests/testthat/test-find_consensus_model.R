context("find_consensus_model")


# Set-up
set.seed(314)
data <- observations_to_histogram(rnorm(10000, mean = 20, sd = 10))
data <- data$histogram_data
models <- fit_distributions(data)

test_that("basic use", {
  
  results <- find_consensus_model(models)
  
  # Formatting
  expect_named(
    results,
    c("jaccard", "intersection", "ks", "mse", "chisq", "consensus")
  )
  
  # Takes the jaccard value in this case
  expect_equal(
    results$jaccard$value,
    results$consensus$value
  )

})

test_that("edge case: reprioritization of metrics", {
  
  results <- find_consensus_model(models, metric = c("intersection", "jaccard", "ks", "mse", "chisq"))
  
  expect_equal(
    results$intersection$value,
    results$consensus$value
  )
  
})

test_that("edge case: one metric", {
  
  models_jaccard <- fit_distributions(data, metric = c("jaccard"))
  results <- find_consensus_model(models_jaccard, metric = c("jaccard"))
  
  expect_named(
    results,
    c("jaccard", "consensus")
  )
  
  expect_equal(
    results$jaccard$value,
    results$consensus$value
  )
  

})

test_that("edge case: one distribution", {
  
  models_norm <- fit_distributions(data, distributions = "norm")
  results <- find_consensus_model(models_norm)
  
  expect_named(
    results,
    c("jaccard", "intersection", "ks", "mse", "chisq", "consensus")
  )
  
  expect_equal(
    results$jaccard$value,
    results$consensus$value
  )
  
  distributions <- sapply(results, `[[`, "dist")
  expect_true(
    all(distributions == "norm")
  )
  
})

test_that("selecting a subset of metrics", {
  
  results <- find_consensus_model(models, metric = c("jaccard", "intersection", "mse"))
  
  expect_named(
    results,
    c("jaccard", "intersection", "mse", "consensus")
  )
  
})

test_that("majority voting - weights", {
  
  # NOTE: FIX THIS IF INVERT METRIC VALUES FOR JACCARD AND INTERSECTION
  models_test <- list(
    "A" = list("metric" = "jaccard", "dist" = "norm", "value" = 0.9),
    "B" = list("metric" = "jaccard", "dist" = "unif", "value" = 0.95),
    "C" = list("metric" = "intersection", "dist" = "norm", "value" = 0.95),
    "D" = list("metric" = "intersection", "dist" = "unif", "value" = 0.9)
  )
  
  results <- find_consensus_model(
    models_test, 
    method = "weighted_majority_vote", 
    metric = c("intersection", "jaccard"),
    weight = c(2, 1)
  )
  
  expect_equal(
    results$intersection$value,
    results$consensus$value
  )
  
  expect_equal(
    results$consensus$dist,
    "norm"
  )
  
})

test_that("rra aggregation", {
  
  # NOTE: RRA doesn't work when there's only 2 distributions - might want to 
  # figure out why and document
  
  models_test <- list(
    "A" = list("metric" = "jaccard", "dist" = "norm", "value" = 0.9),
    "B" = list("metric" = "jaccard", "dist" = "unif", "value" = 0.95),
    "C" = list("metric" = "intersection", "dist" = "norm", "value" = 0.95),
    "D" = list("metric" = "intersection", "dist" = "unif", "value" = 0.9)
  )
  
  expect_error(
    find_consensus_model(models_test, method = "rra", metric = c("intersection", "jaccard")),
    "Ties exist between distributions chosen by metric."
  )
  
  
  models_test <- list(
    "A" = list("metric" = "jaccard", "dist" = "norm", "value" = 0.9),
    "B" = list("metric" = "jaccard", "dist" = "unif", "value" = 0.99),
    "C" = list("metric" = "jaccard", "dist" = "gamma", "value" = 0.999),
    "D" = list("metric" = "mse", "dist" = "unif", "value" = 0.1),
    "E" = list("metric" = "mse", "dist" = "norm", "value" = 0.01),
    "F" = list("metric" = "mse", "dist" = "gamma", "value" = 0.01)
  )
  
  results <- find_consensus_model(
    models_test, 
    method = "rra",
    metric = c("jaccard", "mse")
  )
  
  expect_equal(
    results$jaccard$value,
    results$consensus$value
  )
  
  expect_equal(
    "gamma",
    results$consensus$dist
  )
})