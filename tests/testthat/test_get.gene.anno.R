context("Testing .get.gene.anno")

# Setting Test Parameters
PARAMETERS = list("GTF" = gtf, "GENE" = "ENSGXX")

# .read.gtf
gtf.path = system.file("extdata", package = "DPDE4PM")
gtf = paste0(gtf.path, "/test.gtf")
annotation = DPDE4PM:::.read.gtf(PARAMETERS)

# Testing .get.gene.anno
geneinfo = DPDE4PM:::.get.gene.anno(PARAMETERS = PARAMETERS, ANNOTATION = annotation)

test_that("Output of .get.gene.anno have the appropriate format",{
  
  # results are a data.frame
  expect_is(geneinfo, "list")
  
  # results have the right names
  expected.names = c("anno", "gene", "chr", "strand", "left", "right", "DNA2RNA", "RNA2DNA", "dna_length", "exome_length")
  expect_identical(names(geneinfo), expected.names)
  
  # annotation exists and is data.frame
  expect_is(geneinfo[['anno']], "data.frame")
  
  # gene is character
  expect_is(geneinfo[['gene']], "character")
  
  # chr is character
  expect_is(geneinfo[["chr"]], "character")
  
  # strand
  expect_is(geneinfo[["strand"]], "character")
  expect_true(geneinfo[["strand"]] %in% c("+", "-", "*"))
  
  # left
  expect_is(geneinfo[["left"]], "numeric")
  
  # right
  expect_is(geneinfo[["right"]], "numeric")
  
  # dna_length
  expect_is(geneinfo[["dna_length"]], "numeric")
  
  # exome_length
  expect_is(geneinfo[["exome_length"]], "numeric")
  
  # DNA2RNA
  expect_is(geneinfo[["DNA2RNA"]], "numeric")
  expect_true(length(geneinfo[["DNA2RNA"]]) == dna_length)
  
  # RNA2DNA
  expect_is(geneinfo[["RNA2DNA"]], "numeric")
  expect_true(length(geneinfo[["RNA2DNA"]]) == exome_length)
  
})