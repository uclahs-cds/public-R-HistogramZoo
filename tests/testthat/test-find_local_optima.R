context("find_local_optima")

test_that("base case", {
  x <- c(1,2,1,4,1)
  x_optima <- find_local_optima(x)
  expect_equal(x_optima$min_ind, c(1,3,5))
  expect_equal(x_optima$max_ind, c(2,4))
})

test_that("edge case: entirely flat", {
  x <- rep(1, 5)

  # Returns midpoint
  x_optima <- find_local_optima(x, flat_endpoints = F)
  expect_equal(x_optima$min_ind, c(3))
  expect_equal(x_optima$max_ind, NULL)

  # Returns endpoint
  x_optima <- find_local_optima(x, flat_endpoints = T)
  expect_equal(x_optima$min_ind, c(1, 5))
  expect_equal(x_optima$max_ind, NULL)
})

test_that("edge case: no optima, staircase function", {
  x <- c(1, 1, 1, 2, 2, 3, 3, 3)

  # Returns midpoint
  x_optima <- find_local_optima(x, flat_endpoints = F)
  expect_equal(x_optima$min_ind, c(2))
  expect_equal(x_optima$max_ind, c(7))

  # Returns endpoint
  x_optima <- find_local_optima(x, flat_endpoints = T)
  expect_equal(x_optima$min_ind, c(1, 3))
  expect_equal(x_optima$max_ind, c(6, 8))

})

test_that("flat elements", {
  x <- c(rep(1, 3), rep(2, 3), rep(3, 3), rep(2, 3), rep(1, 3))

  # Returns midpoint
  x_optima <- find_local_optima(x, flat_endpoints = F)
  expect_equal(x_optima$min_ind, c(2, 14))
  expect_equal(x_optima$max_ind, c(8))

  # Returns endpoint
  x_optima <- find_local_optima(x, flat_endpoints = T)
  expect_equal(x_optima$min_ind, c(1, 3, 13, 15))
  expect_equal(x_optima$max_ind, c(7, 9))
})

test_that("threshold works for basic examples", {
  x <- c(1, 2, 1, 2, 1, 3, 5, 3, 1)
  x_optima <- find_local_optima(x, threshold = 1)
  expect_equal(x_optima$min_ind, c(1, 9))
  expect_equal(x_optima$max_ind, c(7))
})

test_that("threshold works for clustered optima", {
  x <- c(1, 4, 1, 2, 1, 4, 1)
  x_optima <- find_local_optima(x, threshold = 1)
  expect_equal(x_optima$min_ind, c(1, 5, 7))
  expect_equal(x_optima$max_ind, c(2, 6))
})

test_that("the one case that threshold doesn't work - endpoints", {
  x <- c(1, 2, 3, 2, 1, 3, 5, 3, 1, 2)
  x_optima <- find_local_optima(x, threshold = 1)
  expect_equal(x_optima$min_ind, c(1))
  expect_equal(x_optima$max_ind, c(7, 10))
})
