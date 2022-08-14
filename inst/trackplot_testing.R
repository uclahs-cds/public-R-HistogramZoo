
library(HistogramZoo)

# Preamble ----------------------------------------------------------------
# Testing trackplot

# Set-up ------------------------------------------------------------------

x = Histogram(c(0, 0, 1, 2, 3, 2, 1, 2, 3, 4, 5, 3, 1, 0))

results = segment_and_fit(x, eps = 0.005)

results_table = summarize_results(results)

# Testing -----------------------------------------------------------------

# Test that you can generate a pretty plot with continuous value column
create_trackplot(
  track_data = results_table,
  row_id = "peak_id",
  metric_id = "value"
)

# Test that you can generate a pretty plot with discrete value column
create_trackplot(
  track_data = results_table,
  row_id = "peak_id",
  metric_id = "peak_id"
)

# Test having 1 row in the results table doesn't screw up orientation
# THIS THROWS AN ERROR - figure out why
# Error in valid.viewport(x, y, width, height, just, gp, clip, mask, xscale,  :
# invalid 'xscale' in viewport
create_trackplot(
  track_data = results_table[1,],
  row_id = "peak_id",
  metric_id = "value",
  xlimits = c(2, 7)
)

# Test that row_id can be character or factor or numeric

# character
create_trackplot(
  track_data = results_table,
  row_id = "dist",
  metric_id = "value"
)

# factor reorders the rows
factor_row_id = results_table
factor_row_id[,"dist"] = factor(factor_row_id[,"dist"], levels = c("norm", "gamma"))
create_trackplot(
  track_data = factor_row_id,
  row_id = "dist",
  metric_id = "value"
)

factor_row_id[,"dist"] = factor(factor_row_id[,"dist"], levels = c("gamma", "norm"))
create_trackplot(
  track_data = factor_row_id,
  row_id = "dist",
  metric_id = "value"
)

# Test that rows can be concatenated
# Weird random error: https://github.com/tidyverse/ggplot2/issues/2514
tmp = create_trackplot(
  track_data = results_table,
  row_id = "region_id",
  metric_id = "value"
)

# Test intervals that don't match xlimits
create_trackplot(
  track_data = results_table,
  row_id = "peak_id",
  metric_id = "value",
  xlimits = c(5, 13)
)

# Test incorrect input
# Incorrect column types
# Incorrect track_data
