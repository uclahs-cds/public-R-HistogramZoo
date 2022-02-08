

library(extraDistr)

# Preamble ----------------------------------------------------------------
# Testing the workflow of the tool on RNA based BED files

# All example data should be in inst/extdata when package is installable
# setwd("/cluster/home/helenzhu/Cluster_Helen/Snakemake_ConsensusPeaks/TestData/dna_bedfiles")
# setwd("inst/extdata")

# Setting up baseline parameters
filenames = list.files(".", pattern = ".bed")
n_fields = 6
gene.or.transcript = "gene"
genes = c("ENSG00000103035.11", "ENSG00000198900.6")

# ### STEP 1 ###
# Loading data into histograms with an appropriate gene model
list.hist = bed.to.hist(
  filenames = filenames,
  n_fields = 6,
  gtf.file = NULL,
  gene.or.transcript = "gene",
  histogram.bin.size = 1
)

# ### STEP 2 ###
# Running the FTC algorithm
ftc.res = bulk.segment.fit(
  coverage.model.obj = list.hist,
  eps = 1,
  seed = NULL,
  truncated.models = TRUE,
  uniform.peak.threshold = 0.75,
  uniform.peak.stepsize = 5,
  remove.low.entropy = T,
  max.uniform = F,
  histogram.metric = c("jaccard", "intersection", "ks", "mse", "chisq")
)

# ### STEP 3 ###
# Compiling results
res = summarize.results(
  coverage.model.obj = ftc.res
)
res

# ### STEP 4 (Optional) ###
# Plotting the results
plot.fitted.segments(
  histogram.names = genes,
  coverage.model.obj = ftc.res,
  file.name = "~/figures/test.pdf"
)

# ### STEP 5 (Optional) ###
# Merging p-values
