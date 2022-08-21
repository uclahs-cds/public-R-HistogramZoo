context("find_consensus_model")


# Set-up
set.seed(314)
data <- observations_to_histogram(rnorm(10000, mean = 20, sd = 10))
data <- data$histogram_data

test_that("base case", {
  
  # Base case
  models <- fit_distributions(data)
  results <- find_consensus_model(models)
  
  # When there's only one metric
  
  # When there's only one dist
  
  
  # Selecting a subset of metrics
  
  
  # Majority voting method
  
  
  # Weighted scores
  
  
  # RRA aggregation method

})