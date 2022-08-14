context("plotting")


# Basic example
x = Histogram(c(0, 0, 1, 2, 3, 2, 1, 2, 3, 4, 5, 3, 1, 0))
results = segment_and_fit(x, eps = 0.005)
results_table = summarize_results(results)

# Long histogram
set.seed(314)
long_x = observations_to_histogram(
  round(
    rnorm(1000, mean = 100, sd = 50)
  ), 
  histogram_bin_width = 5
)
results_long = segment_and_fit(long_x, eps = 1)

# Genomic histogram
genomic_hist = GenomicHistogram(
  c(0, 0, 1, 2, 3, 2, 1, 2, 3, 4, 5, 3, 1, 0), 
  chr = "chr1", 
  strand = c("-"))
genomic_results = segment_and_fit(genomic_hist, eps = 0.005)

# Interval start != interval end
interval_test_data = Histogram(
  histogram_data = c(0, 0, 1, 2, 3, 2, 1, 2, 3, 4, 5, 3, 1, 0),
  interval_start = seq(0, 26, 2),
  interval_end = seq(1, 27, 2))
interval_test_results = segment_and_fit(interval_test_data, eps = 0.005)

# Not all distributions are used
results_dist_subset = segment_and_fit(x, eps = 0.005, distributions = c("unif"))

test_that("create_coverageplot works ", {
  
  # 0. No fitted data
  expect_error(
    create_coverageplot(x),
    NA
  )

  # 1. Basic example
  expect_error(
    create_coverageplot(results),
    NA
  )
  
  # 2. Long example
  expect_error(
    create_coverageplot(results_long),
    NA
  )
  
  # 3. Genomic example
  expect_error(
    create_coverageplot(genomic_results),
    NA
  )
  
  # 4. Interval start !+ interval end example
  expect_error(
    create_coverageplot(interval_test_results),
    NA
  )
  
  # 5. Not all distributions are used
  expect_error(
    create_coverageplot(results_dist_subset),
    NA
  )
  
  # 6. Plot different metric than 'consensus'
  expect_error(
    create_coverageplot(results, "jaccard"),
    NA
  )
  
  # 7. Select a random BPG parameter and modify it
  expect_error(
    create_coverageplot(results, xlab.cex = 1),
    NA
  )
  
})

test_that("create_residualplot works ", {
  
  # 1. Basic example
  expect_error(
    create_residualplot(results),
    NA
  )
  
  # 2. Long example
  expect_error(
    create_residualplot(results_long),
    NA
  )
  
  # 3. Genomic example
  expect_error(
    create_residualplot(genomic_results),
    NA
  )
  
  # 4. Interval start !+ interval end example
  expect_error(
    create_residualplot(interval_test_results),
    NA
  )
  
  # 5. Not all distributions are used
  expect_error(
    create_residualplot(results_dist_subset),
    NA
  )
  
  # 6. Plot different metric than 'consensus'
  expect_error(
    create_residualplot(results, "jaccard"),
    NA
  )
  
  # 7. Select a random BPG parameter and modify it
  expect_error(
    create_residualplot(results, xlab.cex = 1),
    NA
  )
  
  # 8. Add changepoint lines
  expect_error(
    create_residualplot(results, add_changepoint_lines = T),
    NA
  )
  
  # 9. Add changepoint lines + other lines
  expect_error(
    create_residualplot(results, add_changepoint_lines = T, abline.v = c(1), abline.h = c(0, 3)),
    NA
  )
  
})

test_that("create_trackplot works ", {
  
  # 1. Continuous `value`
  expect_error(
    create_trackplot(
      track_data = results_table,
      row_id = "peak_id",
      metric_id = "value"
    ),
    NA
  )
  
  # 2. Discrete `value`
  expect_error(
    create_trackplot(
      track_data = results_table,
      row_id = "peak_id",
      metric_id = "peak_id"
    ),
    NA
  )
  
  # 3. `row` is factor
  factor_row_id = results_table
  factor_row_id[,"dist"] = factor(factor_row_id[,"dist"], levels = c("norm", "gamma"))
  
  expect_error(
    create_trackplot(
      track_data = factor_row_id,
      row_id = "dist",
      metric_id = "value"
    ),
    NA
  )
  
  # 4. `row` is character
  expect_error(
    create_trackplot(
      track_data = results_table,
      row_id = "dist",
      metric_id = "value"
    ),
    NA
  )
  
  # 5. Concatenate rows in tracks
  expect_error(
    create_trackplot(
      track_data = results_table,
      row_id = "region_id",
      metric_id = "value"
    ),
    NA
  )
  # 6. Intervals exceeding x_limits
  expect_warning(
    create_trackplot(
      track_data = results_table,
      row_id = "peak_id",
      metric_id = "value",
      xlimits = c(5, 13)
    ),
    "some intervals are truncated by xlimits"
  )
  
})


test_that("create_layerplot works ", {
  
  cvg_plt = create_coverageplot(results_long, legend = NULL)
  resid_plt = create_residualplot(results_long)
  
  # 1. Testing that this works
  expect_error(
    create_layerplot(
      plot.objects = list(cvg_plt, resid_plt)
    ),
    NA
  )
  
  # 2. Testing that other parameters work
  expect_error(
    create_layerplot(
      plot.objects = list(cvg_plt, resid_plt),
      xlab.label = "TEST",
      xlab.cex = 1
    ),
    NA
  )
  
  
})