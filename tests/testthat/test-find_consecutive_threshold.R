context("find_consecutive_threshold")

test_that("find_consecutive_threshold works", {
  test1 <- c(0,0,0,1,1,1,0,0,0,1,1,1,0,0)
  test2 <- c(1,1,1)

  test1_results <- find_consecutive_threshold(test1)
  test2_results <- find_consecutive_threshold(test2)

  expect_equal(test1_results$start, c(4,10))
  expect_equal(test1_results$end, c(6,12))

  expect_equal(test2_results$start, 1)
  expect_equal(test2_results$end, 3)
})
