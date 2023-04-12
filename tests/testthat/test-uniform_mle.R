test_that("uniform_mle works", {

  expected_zero_mle <- c(0, 1, 50)

  for(x in expected_zero_mle){
    expect_equal(
      uniform_mle(x, x.start = 0, x.end = 1, a = 0, b = 1, log = T),
      0
      )
    }

  # Uniform MLE works for histograms
  x <- Histogram(c(2,2,3,3,4))
  a <- -2
  b <- 10
  res <- uniform_mle(x$histogram_data, x.start = 0, x.end = 5, a = -2, b = 10)

  shift <- 5
  x_shift <- x
  x_shift$interval_start <- x_shift$interval_start + shift
  x_shift$interval_end <- x_shift$interval_end + shift
  # Constant shift should not change results
  expect_equal(res, uniform_mle(x_shift$histogram_data, x.start = head(x_shift$interval_start, 1), x.end = tail(x_shift$interval_end, 1), a = a + shift, b = b + shift))

  # The maximum likelihood should be [min, max]
  expected_min <- head(x$interval_start, n = 1)
  expected_max <- tail(x$interval_end, n = 1)

  a.range <- seq(expected_min - 5, expected_min + 5, by = 1)
  b.range <- seq(expected_max - 5, expected_max + 5, by = 1)
  param.space <- expand.grid(a = a.range, b = b.range)
  neg.log.lik <- mapply(
    FUN = function(a, b) - uniform_mle(x$histogram_data, x.start = head(x$interval_start, 1), x.end = tail(x$interval_end, 1), a, b),
    a = param.space$a,
    b = param.space$b
    )
  expect_equal(
    # The minimum -log-likelihood should be the min and max
    as.numeric(param.space[which.min(neg.log.lik), ]),
    c(expected_min, expected_max)
  )
})
