test_that("find_local_optima works on non consecutive numbers", {
  x <- c(1,2,1,4,1)
  x_optima <- find_local_optima(x)
  expect_equal(x_optima$min_ind, c(1,3,5))
  expect_equal(x_optima$max_ind, c(2,4))
})


test_that("find_local_optima works on consecutive starting numbers", {
  x <- c(1, 1,2,1,4,1)
  x_optima <- find_local_optima(x)
  expect_equal(x_optima$min_ind, c(2,4,6))
  expect_equal(x_optima$max_ind, c(3,5))

  x <- c(1,1,1, 1,2,1,4,1)
  x_optima <- find_local_optima(x)
  expect_equal(x_optima$min_ind, c(4,6,8))
  expect_equal(x_optima$max_ind, c(5,7))
})

test_that("find_local_optima works on consecutive ending numbers", {
  x <- c(1,2,1,4,1, 1)
  x_optima <- find_local_optima(x)
  expect_equal(x_optima$min_ind, c(1,3,5))
  expect_equal(x_optima$max_ind, c(2,4))

  x <- c(1,2,1,4,1, 1,1,1)
  x_optima <- find_local_optima(x)
  expect_equal(x_optima$min_ind, c(1,3,5))
  expect_equal(x_optima$max_ind, c(2,4))
})

test_that("find_local_optima works on consecutive starting and ending numbers", {
  x <- c(1, 1,2,1,4,1, 1)
  x_optima <- find_local_optima(x)
  expect_equal(x_optima$min_ind, c(2,4,6))
  expect_equal(x_optima$max_ind, c(3,5))

  x <- c(1,1,1,1, 1,2,1,4,1, 1,1,1)
  x_optima <- find_local_optima(x)
  expect_equal(x_optima$min_ind, c(5,7,9))
  expect_equal(x_optima$max_ind, c(6,8))
})

test_that("find_local_optima works on even consecutive numbers in middle of sequence", {
  x <- c(1,2,2,1,1,4,4, 1)
  x_optima <- find_local_optima(x)
  expect_equal(x_optima$min_ind, c(1,4,8))
  expect_equal(x_optima$max_ind, c(2,6))
})


# NOTE: This may not be the behavior that we want but is one way of dealing with the consecutive sequences
test_that("find_local_optima works on odd consecutive numbers in middle of sequence", {
  x <- c(1,2,2,2,1,1,1,4,4,4, 1)
  x_optima <- find_local_optima(x)
  expect_equal(x_optima$min_ind, c(1,5,11))
  expect_equal(x_optima$max_ind, c(2,8))

  x <- c(4,2,2,2,1,1,1,4,4,4,5)
  x_optima <- find_local_optima(x)
  expect_equal(x_optima$min_ind, c(2,5, 10))
  expect_equal(x_optima$max_ind, c(1,4,8,11))
})
