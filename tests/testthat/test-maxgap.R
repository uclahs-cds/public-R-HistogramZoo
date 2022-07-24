test_that("find_all_meaningful_gap works on ", {
  x = c(0L, 0L, 1L, 2L, 6L, 11L, 3L, 2L, 0L, 0L, 0L, 2L)
  change.points = c(1, 2, 3, 4, 5, 7, 8, 9, 11, 12)

  meaningful.gaps = find_all_meaningful_gap(x, change.points)

  expect_equal(meaningful.gaps$start, c(1, 9))
  expect_equal(meaningful.gaps$end, c(2, 11))
})
