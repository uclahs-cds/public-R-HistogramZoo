
library(ConsensusPeaks)


# Histogram Example -------------------------------------------------------

x = Histogram(c(0, 0, 1, 2, 3, 2, 1, 2, 3, 4, 5, 3, 1, 0))

results = segment_and_fit(x, eps = 0.005)

results_table = summarize_results(results)

create_coverageplot(
  x
)

create_coverageplot(
  results
)

create_coverageplot(
  results,
  xlab.label = "TEST",
  col = rainbow(4)
)

create_residualplot(
  results
)

# GenomicHistogram Example ------------------------------------------------


x = GenomicHistogram(
  c(0, 0, 1, 2, 3, 2, 1, 2, 3, 4, 5, 3, 1, 0), 
  chr = "chr1", 
  strand = c("-"))

results = segment_and_fit(x, eps = 0.005)

results_table = summarize_results(results)

create_coverageplot(
  results
)

x = GenomicHistogram(
  runif(12), 
  interval_start = seq(2, 120, 10), 
  interval_end = seq(5, 120, 10), 
  chr = "chr1", 
  strand = "-")
