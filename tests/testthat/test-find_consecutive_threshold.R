test_that("find_consecutive_threshold works", {
  test1 <- c(0,0,0,1,1,1,0,0,0,1,1,1,0,0)
  test2 <- c(1,1,1)

  test1.results <- find_consecutive_threshold(test1)
  test2.results <- find_consecutive_threshold(test2)

  expect_equal(test1.results$start, c(4,10))
  expect_equal(test1.results$end, c(6,12))

  expect_equal(test2.results$start, 1)
  expect_equal(test2.results$end, 3)
})
