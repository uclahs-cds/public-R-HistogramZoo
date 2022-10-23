
library(HistogramZoo)

set.seed(314)
my_data <- c(
  rnorm(8000, mean = 50, sd = 2),
  rgamma(10000, shape = 2, rate = 0.4),
  rgamma_flip(9000, shape = 2, rate = 0.4) + 80
)

my_histogram <- observations_to_histogram(my_data)

# A basic coverageplot to visualize the data
create_coverageplot(
  my_histogram,
  main = "Histogram",
  main.cex = 1, 
  xaxis.tck = 0.5, 
  yaxis.tck = 0.5, 
  xlab.cex = 1, 
  ylab.cex = 1, 
  xaxis.cex = 0.8, 
  yaxis.cex = 0.8
)


# Conducting analysis
histogram_fit <- segment_and_fit(
  histogram_obj = my_histogram,
  eps = 1,
  seed = 314,
  truncated_models = TRUE
)

# summarize_results to generate
res <- summarize_results(
  histogram_fit
)

# Histogram Fit
coverage_plt <- create_coverageplot(
  histogram_fit,
  main = "HistogramFit",
  main.cex = 1, 
  xaxis.tck = 0.5, 
  yaxis.tck = 0.5, 
  xlab.cex = 1, 
  ylab.cex = 1, 
  xaxis.cex = 0.8, 
  yaxis.cex = 0.8,
  legend = NULL
)

residual_plt <- create_residualplot(
  histogram_fit,
  main.cex = 0, 
  xaxis.tck = 0.5, 
  yaxis.tck = 0.5, 
  xlab.cex = 0, 
  ylab.cex = 1, 
  xaxis.cex = 0.8, 
  yaxis.cex = 0.8
)

track_plt <- create_trackplot(
  res,
  row_id = "region_id",
  metric_id = "value",
  start_id = "histogram_start",
  end_id = "histogram_end",
  xat = generate_xlabels(histogram_fit, return_xat = T),
  xaxis.lab = generate_xlabels(histogram_fit),
  main.cex = 0, 
  xaxis.tck = 0.5, 
  yaxis.tck = 0.5, 
  xlab.cex = 1, 
  ylab.cex = 1, 
  xaxis.cex = 0.8, 
  yaxis.cex = 0.8
)

create_layerplot(
  plot.objects = list(coverage_plt, residual_plt, track_plt)
)
