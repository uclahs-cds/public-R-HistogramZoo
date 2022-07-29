context("GTF_to_GRangesList")

test_that("GTF_to_GRangesList yields correct results", {

  gtf = system.file("extdata", "genes.gtf", package = "ConsensusPeaks")

  # Basic input
  gtf.gr = GTF_to_GRangesList(gtf)
  expect_named(
    gtf.gr,
    c("ENSG00000178951.9",
      "ENSG00000185129.7"))

  # Selecting transcript
  gtf.transcript = GTF_to_GRangesList(
    gtf,
    gene_or_transcript = "transcript")
  expect_named(
    gtf.transcript,
    c("ENST00000322357.9",
      "ENST00000331327.5",
      "ENST00000502351.1",
      "ENST00000505703.2",
      "ENST00000601588.1",
      "ENST00000651386.1"))

  # Selecting chromosome
  gtf.chr = GTF_to_GRangesList(
    gtf,
    select_chrs = "chr5")
  expect_named(
    gtf.chr,
    c("ENSG00000185129.7")
  )

  # Selecting strand
  gtf.pos = GTF_to_GRangesList(
    gtf,
    select_strand = "+")
  expect_named(
    gtf.pos,
    c("ENSG00000185129.7")
  )

  gtf.neutral = GTF_to_GRangesList(
    gtf,
    select_strand = "*")
  expect_named(
    gtf.neutral,
    c("ENSG00000178951.9",
      "ENSG00000185129.7"))

  # Selecting gene_id
  gtf.gene.select = GTF_to_GRangesList(
    gtf,
    select_ids = "ENSG00000178951.9")
  expect_named(
    gtf.gene.select,
    c("ENSG00000178951.9"))

  # Selecting transcript_id
  gtf.transcript.select = GTF_to_GRangesList(
    gtf,
    gene_or_transcript = 'transcript',
    select_ids = "ENST00000322357.9")
  expect_named(
    gtf.transcript.select,
    c("ENST00000322357.9"))

})
