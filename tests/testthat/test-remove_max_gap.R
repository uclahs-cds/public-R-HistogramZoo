test_that("remove.max.gap works", {
  # set.seed(26)
  # obs.to.int.hist(c(rnorm(100, mean = 1), runif(50, min = 25, max = 35)))
  x = c(`-3` = 0L, `-2` = 0L, `-1` = 3L, `0` = 18L, `1` = 33L, `2` = 30L,
                 `3` = 12L, `4` = 4L, `5` = 0L, `6` = 0L, `7` = 0L, `8` = 0L,
                 `9` = 0L, `10` = 0L, `11` = 0L, `12` = 0L, `13` = 0L, `14` = 0L,
                 `15` = 0L, `16` = 0L, `17` = 0L, `18` = 0L, `19` = 0L, `20` = 0L,
                 `21` = 0L, `22` = 0L, `23` = 0L, `24` = 0L, `25` = 0L, `26` = 7L,
                 `27` = 5L, `28` = 2L, `29` = 7L, `30` = 7L, `31` = 4L, `32` = 2L,
                 `33` = 9L, `34` = 5L, `35` = 2L, `36` = 0L)

  chgpts = find.stepfunction.chgpts(x)


  p = c(3, 8)
  segs = c(start = 3, end = 8)
  p.init = unname(c(segs['start'], chgpts[chgpts > segs['start'] & chgpts < segs['end']], segs['end']))
  p.init = sort(p.init)
  p = ftc.helen(x, p.init, eps = 1)

  mgaps =  meaningful.gaps.local(x = x, seg.points = p, change.points = p.init, min.gap = 1)
  mgaps
  p.pairs = remove.max.gaps.agnostic(p = p, max.gaps = mgaps, remove.short.segment = 1)
  p.pairs
})

test_that("remove.max.gap works on multiple points", {
  p = c(1, 5, 10)
  mgaps = data.frame("start" = 7, "end" = 10)
  p.pairs = remove.max.gaps.agnostic(p = p, max.gaps = mgaps, remove.short.segment = 1)

  expect_equal(p.pairs$start, c(1,5))
  expect_equal(p.pairs$end, c(4,6))

  # Should also work with an empty data frame
  p.pairs.empty = remove.max.gaps.agnostic(p = p, max.gaps = data.frame(), remove.short.segment = 1)

  expect_equal(p.pairs.empty$start, c(1,5))
  expect_equal(p.pairs.empty$end, c(4,10))

  # Short segments should be removed
  mgaps.short = data.frame("start" = 6, "end" = 8)
  p.pairs.no.remove = remove.max.gaps.agnostic(p = p, max.gaps = mgaps.short, remove.short.segment = 0)
  p.pairs.remove = remove.max.gaps.agnostic(p = p, max.gaps = mgaps.short, remove.short.segment = 1)
  expect_equal(nrow(p.pairs.no.remove), 3)
  expect_equal(nrow(p.pairs.remove), 2)
})
