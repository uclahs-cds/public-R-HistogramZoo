test_that("segment.fit works", {
  set.seed(26)
  x.norm.mix = obs.to.int.hist(c(rnorm(100, mean = 1), runif(50, min = 25, max = 35)))
  res = segment.fit.agnostic(x.norm.mix, remove.low.entropy = F)
  res.trunc = segment.fit.agnostic(x.norm.mix, truncated.models = TRUE, remove.low.entropy = F)

  # Expect 2 peaks
  expect_equal(nrow(res$p), 2)
  expect_equal(nrow(res.trunc$p), 2)
})
