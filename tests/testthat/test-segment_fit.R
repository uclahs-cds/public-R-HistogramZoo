test_that("segment.fit works", {
  set.seed(26)
  x.norm.mix = obs.to.int.hist(c(rnorm(100, mean = 1), runif(50, min = 25, max = 35)))
  res = segment.fit.agnostic(x.norm.mix)
  res.trunc = segment.fit.agnostic(x.norm.mix, truncated.models = TRUE, remove.low.entropy =T)

  # Expect 2 peaks
  expect_equal(sum(res$final), 2)
  expect_equal(sum(res.trunc$final), 2)

  # Expect Normal, Uniform
  dists_expect = c("norm", "unif")
  expect_equal(res$dist[res$final == 1], dists_expect)
  expect_equal(res.trunc$dist[res.trunc$final == 1], dists_expect)
  # Expect segment start/ends
})
