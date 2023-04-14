context("utils")

test_that("weighted moments work with constant weights", {
  x <- rnorm(10)
  weights <- rep(1, length(x))
  # Check that weighted versions equal the standard functions
  expect_equal(sd(x), weighted.sd(x, weights))
  expect_equal(var(x), weighted.var(x, weights))
  # moments::skewness computes with biased sd
  expect_equal(moments::skewness(x), weighted.skewness(x, weights, type = 'g'))
})

.replicate_vals <- function(h, x) {
  unlist(mapply(replicate, n = h, expr = x))
}

test_that("weighted moments work with counts", {
  counts <- c(1,2,3,4,1)
  midpoints <- c(42,48,52,58,72)

  long_values <- .replicate_vals(counts, midpoints)

  expect_equal(sd(long_values), weighted.sd(midpoints, counts))
  expect_equal(var(long_values), weighted.var(midpoints, counts))
  expect_equal(moments::skewness(long_values), weighted.skewness(midpoints, counts, type = 'g'))
})


test_that("weighted moments work with on", {
  set.seed(13)
  x <- rnorm(100, sd = 5)
  hist_x <- observations_to_histogram(x, histogram_bin_width = 2)
  mp <- (hist_x$interval_end + hist_x$interval_start) / 2
  long_values <- .replicate_vals(hist_x$histogram_data, mp)
  expect_equal(mean(long_values), weighted.mean(hist_x))
  expect_equal(var(long_values), weighted.var(hist_x))
  expect_equal(sd(long_values), weighted.sd(hist_x))
  expect_equal(moments::skewness(long_values), weighted.skewness(hist_x, type = 'g'))
})
