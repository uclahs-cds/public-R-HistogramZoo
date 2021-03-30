context("Testing .read.gtf")

# Setting Test Parameters
PARAMETERS = list("GTF" = gtf, "GENE" = "ENSGXX")

# Testing .read.gtf
gtf.path = system.file("extdata", package = "DPDE4PM")
gtf = paste0(gtf.path, "/test.gtf")
annotation = DPDE4PM:::.read.gtf(PARAMETERS)

test_that("Output of .read.gtf have the appropriate format",{

  # results are a data.frame
  expect_is(annotation, "data.frame")

  # results have the right column names
  expected.cols = c("chr", "feature", "start", "stop", "strand", "gene", "transcript")
  expect_identical(colnames(annotation), expected.cols)

  # chromosome is character
  expect_is(annotation[,"chr"], "character")

  # function only extracts exons
  expect_true(all(annotation["feature"] == "exon"))

  # start is numeric
  expect_is(annotation[,"start"], "numeric")

  # end is numeric
  expect_is(annotation[,"end"], "numeric")

  # strand
  expect_is(annotation[,"strand"], "character")
  expect_true(all(annotation["strand"] %in% c("+", "-", "*")))

  # gene
  expect_is(annotation[,"gene"], "character")

  # transcript
  expect_is(annotation[,"transcript"], "character")

})
