context("max_gap")

test_that("meaningful_gaps_local works ", {
  x = c(0L, 0L, 1L, 2L, 6L, 11L, 6L, 4L, 2L, 1L, 3L, 5L, 7L, 3L, 1L, 
        0L, 0L, 0L, 1L, 3L, 2L, 0L, 0L, 0L)
  seg.points = c(3, 15, 19, 21)
  change.points = c(2:5, 7:12, 14:16, 18:19, 21:22)
  min.gap = 1

  meaningful.gaps = meaningful_gaps_local(
    x, 
    seg.points,
    change.points,
    min.gap
  )

  expect_equal(meaningful.gaps$start, 16)
  expect_equal(meaningful.gaps$end, 18)
})

test_that("find_all_meaningful_gap works ", {
  x = c(0L, 0L, 1L, 2L, 6L, 11L, 3L, 2L, 0L, 0L, 0L, 2L)
  change.points = c(1, 2, 3, 4, 5, 7, 8, 9, 11, 12)
  
  meaningful.gaps = find_all_meaningful_gap(x, change.points)
  
  expect_equal(meaningful.gaps$start, c(1, 9))
  expect_equal(meaningful.gaps$end, c(2, 11))
})
