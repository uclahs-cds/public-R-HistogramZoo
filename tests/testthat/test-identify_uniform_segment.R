context("identify_uniform_segment")

test_that("identify_uniform_segment works ", {

  # identify_uniform_segments returns output with the correct format
  x <- c(2L, 3L, 2L, 7L, 6L, 7L, 6L, 6L, 7L, 6L, 8L, 7L)

  res <- identify_uniform_segment(
    x = x,
    metric = "jaccard",
    threshold = 0.8,
    stepsize = 1,
    max_sd_size = 0
  )

  expect_named(res, c("par", "dist", "metric", "value", "dens", "seg_start", "seg_end"))

  # Basic reproducibility
  expect_s3_class(res, "ModelFit")
  expect_equal(res[["seg_start"]], 2)
  expect_equal(res[["seg_end"]], 12)
  expect_true(is.null(res[["par"]]))
  expect_equal(res[['dist']], "unif")
  expect_equal(res[['metric']], "jaccard")
  expect_equal(res[['dens']](3, scale = F), 1/11)

  # Sanity check - can this function return the original segment if it's uniform?
  x <- rep(2L, 10)

  res <- identify_uniform_segment(
    x = x,
    metric = "jaccard",
    threshold = 0.8,
    stepsize = 1,
    max_sd_size = 0
  )

  expect_equal(res[["seg_start"]], 1)
  expect_equal(res[["seg_end"]], 10)

})
