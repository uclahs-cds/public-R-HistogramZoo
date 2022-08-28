context("observations_to_histogram")

test_that("Varying inputs to obervations_to_histogram", {

  dataset <- c(rep(1, 3), rep(2, 2), rep(3, 5), rep(4, 1))
  x <- observations_to_histogram(dataset)

  # Returns an object of Histogram-class
  expect_true(inherits(x, "Histogram"))

  # Correct length of bins
  expect_length(x, 4)

  # Correct counts in bins
  expect_equal(x$histogram_data, c(3, 2, 5, 1))

  # Testing appropriate binning
  set.seed(314)
  dataset <- sample(1:100, 1000, replace = T)
  x <- observations_to_histogram(dataset, histogram_bin_width = 10)
  expect_length(x, 10)

  dataset <- sample(1:100, 1000, replace = T)
  x <- observations_to_histogram(dataset, histogram_bin_width = 5)
  expect_length(x, 20)

})
