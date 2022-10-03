test_that("summarize_results returns an appropriate table for a Histogram", {

  x <- Histogram(c(0,0,0,1,1,2,2,3,4,3,2,2,1,0,0,0,1,1,1,0,0))

  sf_results <- segment_and_fit(
    histogram_obj = x,
    histogram_count_threshold = 0,
    eps = 1,
    seed = NULL,
    truncated_models = FALSE,
    uniform_threshold = 0.75,
    uniform_stepsize = 5,
    remove_low_entropy = FALSE,
    min_gap_size = 2,
    min_segment_size = 2,
    max_uniform = FALSE,
    metric = c("jaccard", "intersection", "ks", "mse", "chisq"),
    distributions = c("norm", "gamma", "unif")
  )

  results <- summarize_results(result = sf_results)

  # Checking the format of the output table
  expect_equal(nrow(results), 2)
  expect_true(all(c("region_id", "segment_id", "start", "end",
    "interval_count", "interval_sizes", "interval_starts", "histogram_start",
    "histogram_end", "value", "metric", "dist") %in% colnames(results)))

  expect_is(results[,"region_id"], "character")
  expect_is(results[,"segment_id"], "integer")
  expect_is(results[,"histogram_start"], "numeric")
  expect_is(results[,"histogram_end"], "numeric")

  # Block intervals
  expect_is(results[,"start"], "integer")
  expect_is(results[,"end"], "integer")
  expect_is(results[,"interval_count"], "integer")
  expect_is(results[,"interval_sizes"], "character")
  expect_is(results[,"interval_starts"], "character")

  # Fitted stats
  expect_is(results[,"value"], "numeric")
  expect_is(results[,"metric"], "character")
  expect_is(results[,"dist"], "character")

})

test_that("summarize_results returns an appropriate table for a GenomicHistogram", {

  x <- GenomicHistogram(
          histogram_data = c(0,0,0,1,1,2,2,3,4,3,2,2,1,0,0,0,1,1,1,0,0),
          chr = "chr1",
          strand="*")

  sf_results <- segment_and_fit(
    histogram_obj = x,
    histogram_count_threshold = 0,
    eps = 1,
    seed = NULL,
    truncated_models = FALSE,
    uniform_threshold = 0.75,
    uniform_stepsize = 5,
    remove_low_entropy = FALSE,
    min_gap_size = 2,
    min_segment_size = 2,
    max_uniform = FALSE,
    metric = c("jaccard", "intersection", "ks", "mse", "chisq"),
    distributions = c("norm", "gamma", "unif")
  )

  results <- summarize_results(result = sf_results)

  # Checking the format of the output table
  expect_equal(nrow(results), 2)
  expect_true(all(c("region_id", "segment_id", "start", "end",
    "interval_count", "interval_sizes", "interval_starts", "histogram_start",
    "histogram_end", "value", "metric", "dist") %in% colnames(results)))

  expect_is(results[,"chr"], "character")
  expect_is(results[,"strand"], "character")

})

test_that("summarize_results is capable returning results for different metrics", {

  x <- GenomicHistogram(
          histogram_data = c(0,0,0,1,1,2,2,3,4,3,2,2,1,0,0,0,1,1,1,0,0),
          chr = "chr1",
          strand="*")

  sf_results <- segment_and_fit(
    histogram_obj = x,
    histogram_count_threshold = 0,
    eps = 1,
    seed = NULL,
    truncated_models = FALSE,
    uniform_threshold = 0.75,
    uniform_stepsize = 5,
    remove_low_entropy = FALSE,
    min_gap_size = 2,
    min_segment_size = 2,
    max_uniform = FALSE,
    metric = c("jaccard", "intersection", "ks", "mse", "chisq"),
    distributions = c("norm", "gamma", "unif")
  )

  results_all <- summarize_results(result = sf_results, model_name = "consensus")
  expect_equal(results_all[1, "metric"], "consensus")

  results_jaccard <- summarize_results(result = sf_results, model_name = "jaccard")
  expect_equal(results_jaccard[1, "metric"], "jaccard")

  results_intersection <- summarize_results(result = sf_results, model_name = "intersection")
  expect_equal(results_intersection[1, "metric"], "intersection")

  results_ks <- summarize_results(result = sf_results, model_name = "ks")
  expect_equal(results_ks[1, "metric"], "ks")

  results_mse <- summarize_results(result = sf_results, model_name = "mse")
  expect_equal(results_mse[1, "metric"], "mse")

  results_chisq <- summarize_results(result = sf_results, model_name = "chisq")
  expect_equal(results_chisq[1, "metric"], "chisq")

})
