
library(HistogramZoo)


# Preamble ----------------------------------------------------------------
# Recently realized that we have a bug with plotting disjoint intervals,
# namely that the histograms look weird. Thus, this is an example script
# that tries to rectify this issue


# Generating Data ---------------------------------------------------------

set.seed(314)
histogram_data <- rnorm(10000, mean = 100, sd = 10)
my_histogram <- observations_to_histogram(histogram_data)

create_coverageplot(my_histogram, type = "h")

# Switching up the coordinates
my_disjoint_histogram <- GenomicHistogram(
  histogram_data = my_histogram$histogram_data,
  # Adding disjoint intervals
  interval_start = c(seq(1, 50, 1), seq(61, 87, 1)), 
  interval_end = c(seq(1, 50, 1), seq(61, 87, 1)),
  # This basically only makes sense in a genomic context LOL
  chr = "chr1",
  strand = "+"
)

create_coverageplot(my_disjoint_histogram, type = "h")

# Segment and Fit ---------------------------------------------------------

res <- segment_and_fit(my_histogram, max_uniform = F, remove_low_entropy = F, eps = 1)
disjoint_res <- segment_and_fit(my_disjoint_histogram, max_uniform = F, remove_low_entropy = F, eps = 1)

# Testing Residual Plot ---------------------------------------------------

create_residualplot(
  res,
  add_changepoint_lines = T,
  abline.h = 0,
  abline.v = NULL,
  abline.col = "grey",
  abline.lwd = 1,
  abline.lty = "dotted"
)

create_residualplot(
  disjoint_res,
  add_changepoint_lines = T,
  abline.h = 0,
  abline.v = NULL,
  abline.col = "grey",
  abline.lwd = 1,
  abline.lty = "dotted"
)

# Testing Fitted Coverage Plot --------------------------------------------

create_coverageplot(
  res
)

create_coverageplot(
  disjoint_res
)