test_that("index.to.start.end works", {
  res1 = index.to.start.end(c(1,5,10))
  res2 = index.to.start.end(c(1,2))
  expect_equal(res1$start, c(1, 6))
  expect_equal(res1$end, c(5, 10))

  expect_equal(res2$start, 1)
  expect_equal(res2$end, 2)

  # Expected errors
  expect_error(index.to.start.end(1))
  expect_error(index.to.start.end(NULL))
})
