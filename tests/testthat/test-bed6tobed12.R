context("bed6tobed12")

test_that("bed6tobed12 yields correct results", {
  
  regions = GenomicRanges::GRanges(
    seqnames = "chr1",
    IRanges::IRanges(start = c(1, 7),
                     end = c(4, 10)),
    strand = "+")
  bed12 = bed6tobed12(regions)
  bed12.df = data.frame(bed12)
  
  # Might change this
  bed12.cols = c(
    "seqnames",
    "start",
    "end",
    "width",
    "strand", 
    "score",
    "thickStart",
    "thickEnd",
    "itemRgb",
    "blockCount",
    "blockSizes",
    "blockStarts"
  )
  
  # Output format
  expect_length(bed12, 1)
  expect_named(bed12.df, bed12.cols)
  
  # Expect error if input is incorrect or weird
  incorrect.strand = GenomicRanges::GRanges(
    seqnames = "chr1",
    IRanges::IRanges(start = c(1, 7),
                     end = c(4, 10)),
    strand = c("+", "-"))
  expect_error(bed6tobed12(incorrect.strand))
  
  incorrect.chrom = GenomicRanges::GRanges(
    seqnames = c("chr1", "chr2"),
    IRanges::IRanges(start = c(1, 7),
                     end = c(4, 10)),
    strand = c("+"))
  expect_error(bed6tobed12(incorrect.chrom))
  
})
