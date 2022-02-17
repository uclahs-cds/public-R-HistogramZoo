context("gtf.to.genemodel")

test_that("gtf.to.genemodel yields correct results", {
  
  gtf.file = system.file("extdata", "genes.gtf", package = "ConsensusPeaks")
  
  # Basic input
  gtf.gr = gtf.to.genemodel(gtf.file)
  expect_named(
    gtf.gr, 
    c("ENSG00000178951.9",
      "ENSG00000185129.7"))
 
  # Selecting transcript
  gtf.transcript = gtf.to.genemodel(
    gtf.file,
    gene.or.transcript = "transcript")
  expect_named(
    gtf.transcript,
    c("ENST00000322357.9",
      "ENST00000331327.5",
      "ENST00000502351.1",
      "ENST00000505703.2",
      "ENST00000601588.1",
      "ENST00000651386.1"))
  
  # Selecting chromosome
  gtf.chr = gtf.to.genemodel(
    gtf.file,
    select.chrs = "chr5")
  expect_named(
    gtf.chr,
    c("ENSG00000185129.7")
  )
  
  # Selecting strand
  gtf.pos = gtf.to.genemodel(
    gtf.file,
    select.strand = "+")
  expect_named(
    gtf.pos,
    c("ENSG00000185129.7")
  )
  
  gtf.neutral = gtf.to.genemodel(
    gtf.file,
    select.strand = ".")
  expect_named(
    gtf.neutral, 
    c("ENSG00000178951.9",
      "ENSG00000185129.7"))
  
  # Selecting gene_id
  gtf.gene.select = gtf.to.genemodel(
    gtf.file,
    select.ids = "ENSG00000178951.9")
  expect_named(
    gtf.gene.select, 
    c("ENSG00000178951.9"))
  
  # Selecting transcript_id
  gtf.transcript.select = gtf.to.genemodel(
    gtf.file,
    gene.or.transcript = 'transcript',
    select.ids = "ENST00000322357.9")
  expect_named(
    gtf.transcript.select, 
    c("ENST00000322357.9"))
  
})
