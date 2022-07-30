context("utils")

test_that("index_to_start_end", {
  
  # Right = TRUE
  res1 = index_to_start_end(c(1,5,10))
  expect_equal(res1$start, c(1, 6))
  expect_equal(res1$end, c(5, 10))
  
  res2 = index_to_start_end(c(1,2))
  expect_equal(res2$start, 1)
  expect_equal(res2$end, 2)

  # Right = FALSE
  res1 = index_to_start_end(c(1,5,10), right = FALSE)
  expect_equal(res1$start, c(1, 5))
  expect_equal(res1$end, c(4, 10))
  
  res2 = index_to_start_end(c(1,2), right = FALSE)
  expect_equal(res2$start, 1)
  expect_equal(res2$end, 2)
  
  # Expected errors
  expect_error(index_to_start_end(1))
  expect_error(index_to_start_end(NULL))
})
