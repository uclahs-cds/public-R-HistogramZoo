

# library(valr)
# library(genomation)
# library(GenomicRanges)

# Preamble ----------------------------------------------------------------
# Testing the workflow of the tool on RNA based BED files

# All example data should be in inst/extdata when package is installable
setwd("/cluster/home/helenzhu/Cluster_Helen/Snakemake_ConsensusPeaks/TestData")

# Setting up baseline parameters
filenames = paste0("rna_bedfiles/Sample.", 1:20, ".bed")
n_fields = 12
gtf.file = "genes.gtf"
gene.or.transcript = "gene"
genes = c("ENSG00000198900.6", "ENSG00000103035.11")

#
eps = 10^-4
seed = NULL
truncated.models = FALSE
uniform.peak.threshold = 0.75
uniform.peak.stepsize = 5
remove.low.entropy = T
max.uniform = T
histogram.metric = c("jaccard", "intersection", "ks", "mse", "chisq")

# ### STEP 1 ###
# Loading data into histograms with an appropriate gene model
list.hist = bed.to.hist(
  filenames = filenames,
  n_fields = 12,
  gtf.file = gtf.file,
  gene.or.transcript = "gene",
  histogram.bin.size = 1
)

names(list.hist)
# [1] "histogram.coverage" "gene.model"         "histogram.bin.size"

# Histogram coverage
head(list.hist$histogram.coverage)

# Gene model
list.hist$gene.model

# Bin size
list.hist$histogram.bin.size
# [1] 1

# ### STEP 2 ###
# Running the FTC algorithm
ftc.res = bulk.segment.fit(
  coverage.model.obj = list.hist,
  eps = 10^-4,
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
res = format.results(
  coverage.model.obj = ftc.res
)

res

# ### STEP 4 (Optional) ###
# Plotting the results
plot.fitted.segments(
  histogram.names = genes,
  coverage.model.obj = ftc.res,
  filename = "~/figures/test.pdf"
)

# ### STEP 5 (Optional) ###
# Merging p-values