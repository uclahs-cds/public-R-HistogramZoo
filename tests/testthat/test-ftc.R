context("ftc")

test_that("ftc reproducibility", {

  # FTC works on expected examples
  x <- c(0L, 0L, 1L, 2L, 6L, 11L, 3L, 2L, 0L, 0L, 0L, 2L)
  change_points <- c(1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12)

  # High EPS
  res <- ftc(x, change_points, eps = 1)
  expect_equal(res, c(1, 12))

  # Low EPS
  res <- ftc(x, change_points, eps = 0.33)
  expect_equal(res, c(1 , 2, 3, 4, 7, 8, 9, 10, 11, 12))

  # More complex example
  # TODO: Do we need a 'complex' test case?
  # Any slight updates to the algorithm will cause this to fail...
  # set.seed(314)
  # x <- rnorm(1000, mean = 10, sd = 5)
  # x <- observations_to_histogram(x)
  # x <- x$histogram_data
  # res <- ftc(x, NULL, eps = 0.33)
  # expect_equal(res, c(1, 2, 3, 4, 5, 31, 32, 33, 34))

})

test_that("ftc: density vs. counts", {

  # density
  x <- c(0L, 0L, 1L, 2L, 6L, 11L, 3L, 2L, 0L, 0L, 0L, 2L)
  x <- x/sum(x)
  change_points <- c(1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12)

  # High EPS
  res <- ftc(x, change_points, eps = 1)
  expect_equal(res, c(1, 12))

  # Low EPS - Note that this doesn't just give the same answer above
  res <- ftc(x, change_points, eps = 0.17)
  expect_equal(res, c(1 , 2, 3, 4, 5, 7, 8, 9, 10, 11, 12))

})

test_that("ftc: edge cases", {

  # FTC returns start and end points in a flat region
  x <- rep(0, 10)

  res <- ftc(x, s = NULL, eps = 1)

  expect_equal(res, c(1, 10))

})
