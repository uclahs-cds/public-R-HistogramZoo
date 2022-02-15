

library(ConsensusPeaks)

# Preamble ----------------------------------------------------------------
# Testing the workflow of the tool on RNA based BED files

# All example data should be in inst/extdata when package is installable
setwd("/cluster/home/helenzhu/Cluster_Helen/Snakemake_ConsensusPeaks/TestData")
# setwd("inst/extdata")

# Setting up baseline parameters
filenames = paste0("rna_bedfiles/Sample.", 1:20, ".bed")
n_fields = 12
gtf.file = "genes.gtf"
gene.or.transcript = "gene"
genes = c("ENSG00000178951.9", "ENSG00000185129.7")

# ### STEP 1 ###
# Loading data into histograms with an appropriate gene model
list.hist = transcript.bed.to.histogram(
  filenames = filenames,
  n_fields = 12,
  gtf.file = gtf.file,
  gene.or.transcript = "gene",
  histogram.bin.size = 1
)

# ### STEP 2 ###
# Running the FTC algorithm
gene.results = bulk.segment.fit(
  coverage.model.obj = list.hist,
  eps = 0.005,
  seed = NULL,
  truncated.models = FALSE,
  uniform.peak.threshold = 0.75,
  uniform.peak.stepsize = 5,
  remove.low.entropy = T,
  max.uniform = F,
  histogram.metric = c("jaccard", "intersection", "ks", "mse", "chisq")
)

# ### STEP 3 ###
# Compiling results
gene.results.summary = summarize.results.bulk(
  result = gene.results,
  output.format = "bed"
)
gene.results.summary

# ### STEP 4 (Optional) ###
# Plotting the results
plot.fitted.segments(
  histogram.names = genes,
  coverage.model.obj = ftc.res,
  file.name = "~/figures/test.pdf"
)

# ### STEP 5 (Optional) ###
# Merging p-values
