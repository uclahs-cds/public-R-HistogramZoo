context("Testing DPDE4PM")

# .read.gtf
gtf.path = system.file("extdata", package = "DPDE4PM")
gtf = paste0(gtf.path, "/test.gtf")

# Loading Data
data(norm.1.samples.15)

# Results
results = DPDE4PM(
  GENE = "ENSGXX",
  PEAKS = norm.1.samples.15,
  GTF = gtf,
  RESOLUTION = 50,
  DP.ITERATIONS = 100,
  WEIGHT.THRESHOLD = 0.2,
  N.SD = 1,
  OUTPUTDIR =".",
  PLOT.RESULT=F,
  WRITE.OUTPUT=F,
  OUTPUT.TAG="",
  ALPHA.PRIORS = c(1,2),
  SEED = 123
)

test_that("Output of DPDE4PM has the appropriate format",{
  
  # results are a data.frame
  expect_is(results, "data.frame")
  
  # results have the right names
  expected.names = c("chr", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount",
                     "blockSizes", "blockStarts")
  expect_identical(colnames(results)[1:12], expected.names)
  
  # chr
  expect_is(results['chr'], "character")
  
  # start
  expect_is(results['start'], "numeric")
  
  # end
  expect_is(results['end'], "numeric")
  
  # name
  expect_is(results['name'], 'character')
  
  # strand
  expect_is(results['strand'], 'character')
  expect_true(results["strand"] %in% c("+", "-", "*"))
  
})