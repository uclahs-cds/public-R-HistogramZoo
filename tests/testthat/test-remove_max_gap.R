test_that("remove.max.gap works a single segment", {
  p = c(3, 8)
  mgaps = data.frame("start" = c(3,7), "end" = c(4, 8))
  p.pairs = remove.max.gaps.agnostic(p = p, max.gaps = mgaps, remove.short.segment = 1)

  expect_equal(p.pairs$start, 5)
  expect_equal(p.pairs$end, 6)

  # Should also work with an empty data frame
  p.pairs.empty = remove.max.gaps.agnostic(p = p, max.gaps = data.frame(), remove.short.segment = 1)

  expect_equal(p.pairs.empty$start, 3)
  expect_equal(p.pairs.empty$end, 8)
})

test_that("remove.max.gap works on multiple points", {
  p = c(1, 5, 10)
  mgaps = data.frame("start" = 8, "end" = 10)
  p.pairs = remove.max.gaps.agnostic(p = p, max.gaps = mgaps, remove.short.segment = 1)

  expect_equal(p.pairs$start, c(1,6))
  expect_equal(p.pairs$end, c(5,7))

  # Should also work with an empty data frame
  p.pairs.empty = remove.max.gaps.agnostic(p = p, max.gaps = data.frame(), remove.short.segment = 1)

  expect_equal(p.pairs.empty$start, c(1,6))
  expect_equal(p.pairs.empty$end, c(5,10))

  # Short segments should be removed
  mgaps.short = data.frame("start" = 8, "end" = 9)
  p.pairs.no.remove = remove.max.gaps.agnostic(p = p, max.gaps = mgaps.short, remove.short.segment = 0)
  p.pairs.remove = remove.max.gaps.agnostic(p = p, max.gaps = mgaps.short, remove.short.segment = 1)
  expect_equal(nrow(p.pairs.no.remove), 3)
  expect_equal(nrow(p.pairs.remove), 2)
})
