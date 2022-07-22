test_that("index_to_start_end works", {
  res1 = index_to_start_end(c(1,5,10))
  res2 = index_to_start_end(c(1,2))
  expect_equal(res1$start, c(1, 6))
  expect_equal(res1$end, c(5, 10))

  expect_equal(res2$start, 1)
  expect_equal(res2$end, 2)

  # Expected errors
  expect_error(index_to_start_end(1))
  expect_error(index_to_start_end(NULL))
})
