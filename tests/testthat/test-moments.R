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


test_that("weighted moments work with Histograms + GenomicHistograms", {
  
  set.seed(13)
  
  # Histogram
  x <- rnorm(100, sd = 5)
  hist_x <- observations_to_histogram(x, histogram_bin_width = 2)
  mp <- find_midpoint(hist_x)
  long_values <- .replicate_vals(hist_x$histogram_data, mp)
  
  expect_equal(mean(long_values), weighted.mean(hist_x))
  expect_equal(var(long_values), weighted.var(hist_x))
  expect_equal(sd(long_values), weighted.sd(hist_x))
  expect_equal(moments::skewness(long_values), weighted.skewness(hist_x, type = 'g'))
  
  # GenomicHistogram
  genomichistogram_x <- GenomicHistogram(
    histogram_data = hist_x$histogram_data,
    # Testing that the shift in interval start/end doesn't influence the calculation
    interval_start = hist_x$interval_start + 10,
    interval_end = hist_x$interval_end + 9,
    chr = "chr1",
    strand = "*",
    consecutive_start = hist_x$interval_start + 1,
    consecutive_end = hist_x$interval_end
  )
  
  mp <- find_midpoint(genomichistogram_x)
  long_values <- .replicate_vals(genomichistogram_x$histogram_data, mp)
  
  expect_equal(mean(long_values), weighted.mean(genomichistogram_x))
  expect_equal(var(long_values), weighted.var(genomichistogram_x))
  expect_equal(sd(long_values), weighted.sd(genomichistogram_x))
  expect_equal(moments::skewness(long_values), weighted.skewness(genomichistogram_x, type = 'g'))
  
})
