test_that("segment.fit works without removing low entropy regions", {
  # set.seed(26)
  # x.norm.mix = obs.to.int.hist(c(rnorm(100, mean = 1), runif(50, min = 25, max = 35)))
  x.norm.mix = c(`-3` = 0L, `-2` = 0L, `-1` = 3L, `0` = 18L, `1` = 33L, `2` = 30L,
                 `3` = 12L, `4` = 4L, `5` = 0L, `6` = 0L, `7` = 0L, `8` = 0L,
                 `9` = 0L, `10` = 0L, `11` = 0L, `12` = 0L, `13` = 0L, `14` = 0L,
                 `15` = 0L, `16` = 0L, `17` = 0L, `18` = 0L, `19` = 0L, `20` = 0L,
                 `21` = 0L, `22` = 0L, `23` = 0L, `24` = 0L, `25` = 0L, `26` = 7L,
                 `27` = 5L, `28` = 2L, `29` = 7L, `30` = 7L, `31` = 4L, `32` = 2L,
                 `33` = 9L, `34` = 5L, `35` = 2L, `36` = 0L)
  res = segment.fit.agnostic(x.norm.mix, remove.low.entropy = F)
  res.trunc = segment.fit.agnostic(x.norm.mix, truncated.models = TRUE, remove.low.entropy = F)

  # Expect 2 peaks
  # Peak one should be [3,8]
  # Peak two should be [30, 39]

  p.expect = matrix(c(3,8, 30, 39), nrow = 2, ncol = 2, byrow = TRUE)
  segment.points.matrix = as.matrix(res$p[, c('start', 'end')])
  segment.points.trunc.matrix = as.matrix(res.trunc$p[, c('start', 'end')])
  expect_equivalent(segment.points.matrix, p.expect)
  expect_equivalent(segment.points.trunc.matrix, p.expect)
})
