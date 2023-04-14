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
    GenomicHistogram(
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

test_that("GenomicHistogram-class introns", {

  # Valid introns - self-computing introns
  x <- GenomicHistogram(
    histogram_data = rep(1, 2),
    interval_start = c(1, 5),
    interval_end = c(3, 7),
    chr = "chr1",
    strand = "+"
  )
  expect_equal(x$intron_start, 4)
  expect_equal(x$intron_end, 4)

  # Valid introns - self-computing introns with user-specified introns
  x <- GenomicHistogram(
    histogram_data = rep(1, 2),
    interval_start = c(1, 5),
    interval_end = c(3, 8),
    chr = "chr1",
    strand = "+",
    intron_start = c(6),
    intron_end = c(6)
  )
  expect_equal(x$intron_start, c(4, 6))
  expect_equal(x$intron_end, c(4, 6))

  x <- GenomicHistogram(
    histogram_data = rep(1, 2),
    interval_start = c(1, 5),
    interval_end = c(3, 8),
    chr = "chr1",
    strand = "+",
    intron_start = c(4, 6),
    intron_end = c(4, 6)
  )
  expect_equal(x$intron_start, c(4, 6))
  expect_equal(x$intron_end, c(4, 6))

  # Valid introns - intron between bins
  expect_error(
    GenomicHistogram(
      histogram_data = rep(1, 2),
      interval_start = c(1, 5),
      interval_end = c(3, 7),
      chr = "chr1",
      strand = "+",
      intron_start = c(4),
      intron_end = c(4)
    ),
    NA
  )

  # Valid introns - introns within bin
  expect_error(
    GenomicHistogram(
      histogram_data = rep(1, 2),
      interval_start = c(1, 5),
      interval_end = c(4, 7),
      chr = "chr1",
      strand = "+",
      intron_start = c(3),
      intron_end = c(3)
    ),
    NA
  )

  # An odd edge case of an intron which is actually acceptable
  expect_error(
    GenomicHistogram(
      histogram_data = rep(1, 2),
      interval_start = c(1, 5),
      interval_end = c(4, 7),
      chr = "chr1",
      strand = "+",
      intron_start = c(4),
      intron_end = c(4)
    ),
    NA
  )

  # Valid introns - invalid bins (i.e. not of equal length)
  # Middle bin shorter because of intron
  expect_error(
    GenomicHistogram(
      histogram_data = rep(1, 3),
      interval_start = c(1, 4, 8),
      interval_end = c(3, 6, 10),
      chr = "chr1",
      strand = "+",
      intron_start = c(5),
      intron_end = c(5)
    )
  )

  # Invalid introns: overlapping introns
  expect_error(
    GenomicHistogram(
      histogram_data = c(1),
      interval_start = c(1),
      interval_end = c(10),
      chr = "chr1",
      strand = "+",
      intron_start = c(4, 5),
      intron_end = c(6, 7)
    )
  )

  # Invalid introns: intron out of range
  expect_error(
    GenomicHistogram(
      histogram_data = c(1),
      interval_start = c(1),
      interval_end = c(10),
      chr = "chr1",
      strand = "+",
      intron_start = c(8),
      intron_end = c(11)
    )
  )

  # Subsetting introns - checking that subsetting a Histogram accurately
  # subsets introns
  x <- GenomicHistogram(
    histogram_data = rep(1, 2),
    interval_start = c(1, 5),
    interval_end = c(3, 7),
    chr = "chr1",
    strand = "+",
    intron_start = c(2, 6),
    intron_end = c(2, 6)
  )

  x_subset <- x[1]
  expect_equal(x_subset$interval_start, 1)
  expect_equal(x_subset$interval_end, 3)
  expect_equal(x_subset$intron_start, 2)
  expect_equal(x_subset$intron_end, 2)

  # Testing last bin bin_width re-estimation
  x <- GenomicHistogram(
    histogram_data = rep(1, 2),
    interval_start = c(1, 3),
    interval_end = c(1, 5), 
    chr = "chr1",
    strand = "+",
    intron_start = c(4),
    intron_end = c(4)
  )

  x_subset <- x[2]
  expect_equal(x_subset$bin_width, 2)

  # Testing print with introns
  x <- GenomicHistogram(
    histogram_data = rep(1, 3),
    interval_start = c(1,4,6),
    interval_end = c(3,5,7),
    chr = "chr1",
    strand = "+",
    intron_start = c(2),
    intron_end = c(2)
  )

  # Print
  # TODO: Add the expected print output
  expect_invisible(print(x))

})
