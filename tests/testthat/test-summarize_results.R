test_that("summarize.results returns an appropriate table for a basic histogram", {

  x <- c(0,0,0,1,1,2,2,3,4,3,2,2,1,0,0,0,1,1,1,0,0)

  histogram.results <- segment.fit.agnostic(
    x = x,
    histogram.count.threshold = 0,
    eps = 1,
    seed = NULL,
    truncated.models = FALSE,
    uniform.peak.threshold = 0.75,
    uniform.peak.stepsize = 5,
    remove.low.entropy = FALSE,
    min.gap.size = 2,
    min.peak.size = 2,
    max.uniform = FALSE,
    histogram.metric = c("jaccard", "intersection", "ks", "mse", "chisq")
  )

  results = summarize.results(
    segment.fit.agnostic.result = histogram.results,
    output.format = "stats.only")

  # Checking the format of the output table
  expect_equal(nrow(results), 2)
  expect_is(results[,"start"], "numeric")
  expect_is(results[,"end"], "numeric")
  expect_is(results[,"value"], "numeric")
  expect_is(results[,"dist"], "character")
  expect_is(results[,"params"], "character")

})
