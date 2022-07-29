
library(ConsensusPeaks)

# Preamble ----------------------------------------------------------------
# Testing bigWig files

# Loading Data ------------------------------------------------------------

filename = system.file("extdata", "bigwigs",  "S1.bw", package = "ConsensusPeaks")
strand = "."
score.threshold = 1

regions = GenomicRanges::GRanges(
  seqnames = "chr1",
  IRanges::IRanges(start = c(17950, 19350),
          end = c(18000, 19600)),
  strand = ".")

histograms = bigwig.to.histogram(
  filename = filename,
  strand = strand
  score.threshold = score.threshold,
  regions = regions,
  gtf.file = NULL,
  histogram.bin.size = 10)

results = bulk.segment.fit(
  coverage.model.obj = histograms,
  eps = 0.005,
  seed = NULL,
  truncated.models = FALSE,
  uniform.peak.threshold = 0.75,
  uniform.peak.stepsize = 5,
  remove.low.entropy = T,
  max.uniform = F,
  histogram.metric = c("jaccard", "intersection", "ks", "mse", "chisq")
)

gene.results.summary = summarize.results.bulk(
  result = results,
  output.format = "bed"
)
