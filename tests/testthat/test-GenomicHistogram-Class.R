context("GenomicHistogram-class")

test_that("Generating a valid Histogram object", {

  x <- GenomicHistogram(
    histogram_data = c(1, 2, 3, 4, 3, 3, 1),
    chr = "chr1"
  )

  # strand
  expect_equal(x$strand, "*")

  # chromosome
  expect_equal(x$chr, "chr1")

  # region_id
  expect_equal(x$region_id, "chr1:1-7:*")

})

test_that("Generating an invalid GenomicHistogram object", {

  expect_error(
    Histogram(
      runif(10),
      interval_start = 1:10,
      interval_end = 1:10 + 1
    )
  )

  # strand
  expect_error(
    GenomicHistogram(
      histogram_data = c(1, 2, 3, 4, 3, 3, 1),
      strand = "TEST"
    )
  )

})

test_that("GenomicHistogram-class methods", {

  x <- GenomicHistogram(
    histogram_data = c(1, 2, 3, 4, 3, 3, 1),
    chr = "chr1"
  )

  # Print
  # TODO: Add the expected print output
  expect_invisible(print(x))

  # Length
  expect_length(x , 7)

  # Extract Indices
  x_subset <- x[2:5]
  expect_equal(
    x_subset$histogram_data,
    c(2, 3, 4, 3)
  )

  # reassign_ids
  x_reassign <- reassign_region_id(x, "TEST")
  expect_equal(x_reassign$region_id, "TEST")

})
