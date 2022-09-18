context("segment_and_fit")


x_data <- c(`-3` = 0L, `-2` = 0L, `-1` = 3L, `0` = 18L, `1` = 33L, `2` = 30L,
             `3` = 12L, `4` = 4L, `5` = 0L, `6` = 0L, `7` = 0L, `8` = 0L,
             `9` = 0L, `10` = 0L, `11` = 0L, `12` = 0L, `13` = 0L, `14` = 0L,
             `15` = 0L, `16` = 0L, `17` = 0L, `18` = 0L, `19` = 0L, `20` = 0L,
             `21` = 0L, `22` = 0L, `23` = 0L, `24` = 0L, `25` = 0L, `26` = 7L,
             `27` = 5L, `28` = 2L, `29` = 7L, `30` = 7L, `31` = 4L, `32` = 2L,
             `33` = 9L, `34` = 5L, `35` = 2L, `36` = 0L)

x_histogram <- Histogram(
  x_data
)

x_genomic <- GenomicHistogram(
  x_data,
  chr = "chr1",
  strand = "-"
)


test_that("base case: segment_and_fit returns expected output", {

  res <- segment_and_fit(x_histogram)
  res_genomic <- segment_and_fit(x_genomic)

  expect_named(
    res,
    c("histogram_data",
      "interval_start",
      "interval_end",
      "region_id",
      # Results
      "models", "p",
      # Parameters
      "optima_threshold",
      "optima_flat_endpoints",
      "histogram_count_threshold",
      "eps",
      "remove_low_entropy",
      "min_gap_size",
      "min_segment_size",
      "seed",
      "max_uniform",
      "uniform_threshold",
      "uniform_stepsize",
      "uniform_max_sd",
      "truncated_models",
      "metric",
      "distributions",
      "consensus_method",
      "metric_weights"
    )
  )

  expect_named(
    res_genomic,
    c("histogram_data",
      "interval_start",
      "interval_end",
      "region_id",
      "chr",
      "strand",
      # Results
      "models", "p",
      # Parameters
      "optima_threshold",
      "optima_flat_endpoints",
      "histogram_count_threshold",
      "eps",
      "remove_low_entropy",
      "min_gap_size",
      "min_segment_size",
      "seed",
      "max_uniform",
      "uniform_threshold",
      "uniform_stepsize",
      "uniform_max_sd",
      "truncated_models",
      "metric",
      "distributions",
      "consensus_method",
      "metric_weights"
    )
  )

  p_expect <- matrix(c(3,8, 30, 39), nrow = 2, ncol = 2, byrow = TRUE)
  segment_points_matrix <- as.matrix(res$p[, c('start', 'end')])
  expect_equivalent(segment_points_matrix, p_expect)

})

test_that("subfunction: find_local_optima", {

  x_uniform <- Histogram(c(1, rep(2, 10)))

  # Thankfully, this doesn't throw an error
  expect_error(
    segment_and_fit(x_uniform, optima_flat_endpoints = F),
    NA
  )

  # This doesn't change anything except make things faster
  res <- segment_and_fit(x_histogram, optima_threshold = 5, optima_flat_endpoints = F)

  expect_equal(res$optima_threshold, 5)
  expect_true(!res$optima_flat_endpoints)

  p_expect <- matrix(c(3,8, 30, 39), nrow = 2, ncol = 2, byrow = TRUE)
  segment_points_matrix <- as.matrix(res$p[, c('start', 'end')])
  expect_equivalent(segment_points_matrix, p_expect)

})

test_that("subfunction: find_consecutive_threshold", {

  res <- segment_and_fit(x_histogram, histogram_count_threshold = 2)

  expect_equal(res$histogram_count_threshold, 2)
  p_expect <- matrix(c(3,8, 33, 35), nrow = 2, ncol = 2, byrow = TRUE)
  segment_points_matrix <- as.matrix(res$p[, c('start', 'end')])
  expect_equivalent(segment_points_matrix, p_expect)

})

test_that("subfunction: ftc", {

  res <- segment_and_fit(x_histogram, eps = 0.005)

  expect_equal(res$eps, 0.005)
  p_expect <- matrix(c(3, 5, 6, 8, 30, 32), nrow = 3, ncol = 2, byrow = TRUE)
  segment_points_matrix <- as.matrix(res$p[, c('start', 'end')])
  expect_equivalent(segment_points_matrix, p_expect)

})

test_that("subfunction: meaningful_gaps_local and remove_max_gaps", {

  # Potentially re-examine this
  res <- segment_and_fit(x_histogram, remove_low_entropy = F)
  res_summary <- summarize_results(res)

  expect_true(!res$remove_low_entropy)
  expect_equal(res_summary$end[2], 38)

})

test_that("subfunction: identify_uniform_segment", {

  x_uniform <- Histogram(c(1, rep(2, 10)))
  res <- segment_and_fit(x_uniform, max_uniform = T, uniform_threshold = 0.75, uniform_stepsize = 1)
  res_summary <- summarize_results(res)

  expect_true(res$max_uniform)
  expect_equal(res$uniform_threshold, 0.75)
  expect_equal(res$uniform_stepsize, 1)
  expect_equal(res_summary[1, "start"], 2)
  expect_equal(res_summary[1, "dist"], "unif")
  expect_equal(res_summary[1, "value"], 1)

})

test_that("subfunction: fit_distributions", {

  res <- segment_and_fit(x_histogram, distributions = c("norm"))
  res_summary <- summarize_results(res)

  expect_equal(res$distributions, "norm")
  expect_true(all(res_summary$dist == "norm"))

  res2 <- segment_and_fit(x_histogram, distributions = c("norm"), truncated_models = T)
  res2_summary <- summarize_results(res2)

  expect_true(res2$truncated_models)
  expect_true(res2_summary$value[2] > res_summary$value[2])
})

test_that("subfunction: find_consensus_model", {

  res <- segment_and_fit(x_histogram, consensus_method = "rra")
  res_summary <- summarize_results(res)

  expect_equal(res$consensus_method, "rra")
  expect_equal(res_summary[2, "dist"], "unif")

  res <- segment_and_fit(x_histogram, distributions = c("norm", "gamma"), metric = c("intersection", "jaccard"), metric_weights = c(5, 1))
  res_summary <- summarize_results(res)

  expect_equal(res$consensus_method, "weighted_majority_vote")
  expect_equal(res_summary[2, "dist"], "gamma")

})
