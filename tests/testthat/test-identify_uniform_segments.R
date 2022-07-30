context("identify_uniform_segment")

test_that("identify_uniform_segment works ", {
  
  # identify_uniform_segments returns output with the correct format
  x <- c(2L, 3L, 2L, 7L, 6L, 7L, 6L, 6L, 7L, 6L, 8L, 7L)
  
  res <- identify_uniform_segment(
    x = x,
    metric = "jaccard",
    threshold = 0.8,
    stepsize = 1,
    max.sd.size = 0
  )
  
  expect_named(res, c("start", "end", "metric", "length"))
  expect_equal(nrow(res), 1)
  expect_true(res[1, "end"] - res[1, "start"] + 1 == res[1, "length"])
  
  # Basic reproducibility
  expect_equal(res[1, "start"], 2)
  expect_equal(res[1, "end"], 12)
  
  # Sanity check - can this function return the original segment if it's uniform?
  x <- rep(2L, 10)
  
  res <- identify_uniform_segment(
    x = x,
    metric = "jaccard",
    threshold = 0.8,
    stepsize = 1,
    max.sd.size = 0
  )
  
  expect_equal(res[1, "start"], 1)
  expect_equal(res[1, "end"], 10)
  
})