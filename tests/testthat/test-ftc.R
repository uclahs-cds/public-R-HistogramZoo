context("ftc")

test_that("ftc works ", {

  # FTC works on expected examples
  x = c(0L, 0L, 1L, 2L, 6L, 11L, 3L, 2L, 0L, 0L, 0L, 2L)
  change.points = c(1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12)
  
  res = ftc(x, change.points, eps = 1)
  
  expect_equal(res, c(1, 12))
  
  # Erronous input to FTC leads to an error
})