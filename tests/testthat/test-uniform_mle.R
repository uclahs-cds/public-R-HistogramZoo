test_that("uniform.mle works", {

  expected_zero_mle <- c(0, 1, 50)

  for(x in expected_zero_mle){
    expect_equal(
      uniform.mle(x, a = 0, b = 1, log = T),
      0
      )
    }

  # Uniform MLE works for histograms
  x <- Histogram(c(2,2,3,3,4))
  a <- -2
  b <- 10
  res <- uniform.mle(x, a = -2, b = 10)

  shift <- 5
  x_shift <- x
  x_shift$interval_start <- x_shift$interval_start + shift
  x_shift$interval_end <- x_shift$interval_end + shift
  # Constant shift should not change results
  expect_equal(res, uniform.mle(x_shift, a + shift, b + shift))

  # The maximum likelihood should be [min, max]
  expected_min <- head(x$interval_start, n = 1)
  expected_max <- tail(x$interval_end, n = 1)

  a.range <- seq(expected_min - 5, expected_min + 5, by = 1)
  b.range <- seq(expected_max - 5, expected_max + 5, by = 1)
  param.space <- expand.grid(a = a.range, b = b.range)
  neg.log.lik <- mapply(
    FUN = function(a, b) - uniform.mle.Histogram(x, a, b),
    a = param.space$a,
    b = param.space$b
    )
  expect_equal(
    # The minimum -log-likelihood should be the min and max
    as.numeric(param.space[which.min(neg.log.lik), ]),
    c(expected_min, expected_max)
  )
})
