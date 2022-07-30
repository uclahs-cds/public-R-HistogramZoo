context("fit_distributions")

test_that("fit_distributions works ", {
  
  # fit_distributions returns accurately formatted results
  metric <- c("jaccard", "intersection", "ks", "mse", "chisq")
  distributions <- c("norm", "gamma", "unif")
  todo <- expand.grid(metric, distributions, stringsAsFactors = F)
  result_names <- paste0(todo[,1], ".", todo[,2])
  fit_names <- c("par", "dist", "metric", "value", "dens")
  
  histogram_data <- round(rnorm(100, mean = 0, sd = 5))
  histogram_data <- table(histogram_data)
  
  res <- fit_distributions(
    histogram_data,
    metric = metric,
    truncated = F,
    distributions = distributions
  )
  
  expect_length(res, 15)
  expect_named(res, result_names)
  expect_named(res[[1]], fit_names)
  
  # TODO: fit_distributions returns the correct result based on input
  
  # TODO: fit_distributions throws an error when the input isn't correct
})