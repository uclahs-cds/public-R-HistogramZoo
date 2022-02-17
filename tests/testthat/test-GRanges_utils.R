context("Testing input and parameter variations to GRanges utils")

granges = GenomicRanges(
  seqnames = "chr1",
  IRanges::IRanges(
    start = c(3, 5),
    end = c(4, 9)
  ),
  strand = "+"
)

test_that("base0tobase1 works", {
  
  base1 = base0tobase1(granges)
  expect_equal(GenomicRanges::start(granges), c(4, 6))
  expect_equal(GenomicRanges::end(granges), c(4, 9))
  
})

test_that("base1tobase0 works", {
  
  base0 = base1tobase0(granges)
  expect_equal(GenomicRanges::start(granges), c(2, 4))
  expect_equal(GenomicRanges::end(granges), c(4, 9))
  
})

test_that("generate.identifiers works", {
  
  id = generate.identifiers(granges)
  expect_length(id, 2)
  expect_equal(id, c("chr1:3-4", "chr1:5-9"))
  
})