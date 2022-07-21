context("Testing input and parameter variations to GRanges utils")

granges = GenomicRanges::GRanges(
  seqnames = "chr1",
  IRanges::IRanges(
    start = c(3, 5),
    end = c(4, 9)
  ),
  strand = "+"
)

test_that("base0.to.base1 works", {

  base1 = base0.to.base1(granges)
  expect_equal(GenomicRanges::start(base1), c(4, 6))
  expect_equal(GenomicRanges::end(base1), c(4, 9))

})

test_that("base1.to.base0 works", {

  base0 = base1.to.base0(granges)
  expect_equal(GenomicRanges::start(base0), c(2, 4))
  expect_equal(GenomicRanges::end(base0), c(4, 9))

})

test_that("generate.identifiers works", {

  id = generate.identifiers(granges)
  expect_length(id, 2)
  expect_equal(id, c("chr1:3-4:+", "chr1:5-9:+"))

})
