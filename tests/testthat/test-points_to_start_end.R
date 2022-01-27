test_that("points.to.start.end works", {
  res1 = points.to.start.end(c(1,5,10))
  res2 = points.to.start.end(c(1,2))
  expect_equal(res1$start, c(1, 5))
  expect_equal(res1$end, c(4, 10))

  expect_equal(res2$start, 1)
  expect_equal(res2$end, 2)

  # Expected errors
  expect_error(points.to.start.end(1))
  expect_error(points.to.start.end(NULL))
})
