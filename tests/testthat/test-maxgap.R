context("max_gap")

test_that("meaningful_gaps_local: basic reproducibility", {
  
  x <- c(0L, 0L, 1L, 2L, 6L, 11L, 6L, 4L, 2L, 1L, 3L, 5L, 7L, 3L, 1L,
        0L, 0L, 0L, 1L, 3L, 2L, 0L, 0L, 0L)
  seg_points <- c(3, 15, 19, 21)
  change_points <- c(2:5, 7:12, 14:16, 18:19, 21:22)
  min_gap <- 1

  meaningful_gaps <- meaningful_gaps_local(
    x,
    seg_points,
    change_points,
    min_gap
  )

  expect_equal(meaningful_gaps$start, 16)
  expect_equal(meaningful_gaps$end, 18)
  
})

test_that("meaningful_gaps_local: edge case", {
  
  x <- rep(0, 10)
  seg_points <- c(1, 10)
  change_points <- 1:10
  
  meaningful_gaps <- meaningful_gaps_local(
    x,
    seg_points,
    change_points,
    min_gap = 1
  )
  
  expect_equal(nrow(meaningful_gaps), 0)
  
})

test_that("find_all_meaningful_gap: basic reproducibility", {
  
  x <- c(0L, 0L, 1L, 2L, 6L, 11L, 3L, 2L, 0L, 0L, 0L, 2L)
  change_points <- c(1, 2, 3, 4, 5, 7, 8, 9, 11, 12)

  meaningful_gaps <- find_all_meaningful_gap(x, change_points)

  expect_equal(meaningful_gaps$start, c(1, 9))
  expect_equal(meaningful_gaps$end, c(2, 11))
  
})


test_that("find_all_meaningful_gap: edge case", {
  
  x <- rep(0, 10)
  change_points <- 1:10
  
  meaningful_gaps <- find_all_meaningful_gap(x, change_points)
  
  expect_equal(nrow(meaningful_gaps), 0)
})
