context("remove_max_gaps")

test_that("remove_max_gaps works a single segment", {
  start_end_points = index_to_start_end( c(3, 8) )
  mgaps = data.frame("start" = c(3,7), "end" = c(4, 8))
  p_pairs = remove_max_gaps(start_end_points = start_end_points, max_gaps = mgaps, remove_short_segment = 1)

  expect_equal(p_pairs$start, 5)
  expect_equal(p_pairs$end, 6)

  # Should also work with an empty data frame
  p_pairs_empty = remove_max_gaps(start_end_points = start_end_points, max_gaps = data.frame(), remove_short_segment = 1)

  expect_equal(p_pairs_empty$start, 3)
  expect_equal(p_pairs_empty$end, 8)
})

test_that("remove_max_gaps works on multiple points", {
  start_end_points = index_to_start_end( c(1, 5, 10) )
  mgaps = data.frame("start" = 8, "end" = 10)
  p_pairs = remove_max_gaps(start_end_points = start_end_points, max_gaps = mgaps, remove_short_segment = 1)

  expect_equal(p_pairs$start, c(1,6))
  expect_equal(p_pairs$end, c(5,7))

  # Should also work with an empty data frame
  p_pairs_empty = remove_max_gaps(start_end_points = start_end_points, max_gaps = data.frame(), remove_short_segment = 1)

  expect_equal(p_pairs_empty$start, c(1,6))
  expect_equal(p_pairs_empty$end, c(5,10))

  # Short segments should be removed
  mgaps_short = data.frame("start" = 8, "end" = 9)
  p_pairs_no_remove = remove_max_gaps(start_end_points = start_end_points, max_gaps = mgaps_short, remove_short_segment = 0)
  p_pairs_remove = remove_max_gaps(start_end_points = start_end_points, max_gaps = mgaps_short, remove_short_segment = 1)
  expect_equal(nrow(p_pairs_no_remove), 3)
  expect_equal(nrow(p_pairs_remove), 2)
})
