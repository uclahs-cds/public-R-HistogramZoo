context("coverage_to_histogram")

# Different types of regions
one_segment_region <- GenomicRanges::GRanges(
  seqnames = "chr1",
  IRanges::IRanges(
    start = c(1),
    end = c(10)
  ),
  strand = "+"
)

two_segment_region <- GenomicRanges::GRanges(
  seqnames = "chr1",
  IRanges::IRanges(
    start = c(1, 21),
    end = c(10, 30)
  ),
  strand = "+"
)

multi_intron_region <- GenomicRanges::GRanges(
  seqnames = "chr1",
  IRanges::IRanges(
    start = c(1, 6, 9),
    end = c(4, 7, 20)
  ),
  strand = "+"
)

# coverage
coverage <- GenomicRanges::coverage(
  GenomicRanges::GRanges(
    seqnames = "chr1",
    IRanges::IRanges(
      start = c(1, 21),
      end = c(20, 30)
    ),
    strand = "+"
  )
)

uneven_coverage <- GenomicRanges::coverage(
  GenomicRanges::GRanges(
    seqnames = "chr1",
    IRanges::IRanges(
      start = c(1, 21, 21),
      end = c(20, 30, 30)
    ),
    strand = "+"
  )
)

# region id
region_id <- "TEST"

test_that("Expected behavior of coverage_to_histogram", {
  
  
  # Bin size 1 - base case
  base_case <- coverage_to_histogram(
    region = one_segment_region,
    coverage = coverage,
    region_id =  region_id,
    histogram_bin_size = 1
  )
  
  expect_equal(
    base_case$histogram_data,
    rep(1, 10)
  )
  
  expect_equal(
    base_case$interval_start,
    base_case$interval_end
  )
  
  expect_equal(
    base_case$intron_start,
    integer(0)
  )
  
  expect_equal(
    base_case$consecutive_start,
    seq(1, 10)
  )
  
  expect_equal(
    base_case$consecutive_end,
    seq(1, 10)
  )
  
  # Checking that the last bin is indeed a different size if length of object is different
  non_divisible_bin_length <- coverage_to_histogram(
    region = one_segment_region,
    coverage = coverage,
    region_id =  region_id,
    histogram_bin_size = 4
  )
  
  expect_equal(
    non_divisible_bin_length$histogram_data,
    rep(1, 3)
  )
  
  expect_equal(
    non_divisible_bin_length$interval_start,
    c(1, 5, 9)
  )

  expect_equal(
    non_divisible_bin_length$interval_end,
    c(4, 8, 10)
  )
  
  expect_equal(
    non_divisible_bin_length$intron_start,
    integer(0)
  )

  expect_equal(
    non_divisible_bin_length$consecutive_start,
    c(1, 5, 9)
  )
  
  expect_equal(
    non_divisible_bin_length$consecutive_end,
    c(4, 8, 10)
  )
  
  # Checking that bins span introns
  base_intron_case <- coverage_to_histogram(
    region = two_segment_region,
    coverage = coverage,
    region_id =  region_id,
    histogram_bin_size = 4
  )
  
  expect_equal(
    base_intron_case$histogram_data,
    rep(1, 5)
  )
  
  expect_equal(
    base_intron_case$interval_start,
    c(1, 5, 9, 23, 27)
  )
  
  expect_equal(
    base_intron_case$interval_end,
    c(4, 8, 22, 26, 30)
  )
  
  expect_equal(
    base_intron_case$intron_start,
    11
  )
  
  expect_equal(
    base_intron_case$intron_end,
    20
  )
  
  expect_equal(
    base_intron_case$consecutive_start,
    c(1, 5, 9, 13, 17)
  )
  
  expect_equal(
    base_intron_case$consecutive_end,
    c(4, 8, 12, 16, 20)
  )
  
  # Multiple introns per bin
  multi_intron_case <- coverage_to_histogram(
    region = multi_intron_region,
    coverage = coverage,
    region_id =  region_id,
    histogram_bin_size = 4
  )

  expect_equal(
    multi_intron_case$histogram_data,
    rep(1, 5)
  )
  
  expect_equal(
    multi_intron_case$interval_start,
    c(1, 6, 11, 15, 19)
  )
  
  expect_equal(
    multi_intron_case$interval_end,
    c(4, 10, 14, 18, 20)
  )
  
  expect_equal(
    multi_intron_case$intron_start,
    c(5, 8)
  )
  
  expect_equal(
    multi_intron_case$intron_end,
    c(5, 8)
  ) 
  
  expect_equal(
    multi_intron_case$consecutive_start,
    c(1, 5, 9, 13, 17)
  )
  
  expect_equal(
    multi_intron_case$consecutive_end,
    c(4, 8, 12, 16, 18)
  )
  
  # Uneven coverage spanning bins
  uneven_coverage_case <- coverage_to_histogram(
    region = two_segment_region,
    coverage = uneven_coverage,
    region_id =  region_id,
    histogram_bin_size = 4
  )
  
  expect_equal(
    uneven_coverage_case$histogram_data,
    c(1, 1, 1.5, 2, 2)
  )
  
  expect_equal(
    uneven_coverage_case$interval_start,
    c(1, 5, 9, 23, 27)
  )
  
  expect_equal(
    uneven_coverage_case$interval_end,
    c(4, 8, 22, 26, 30)
  )
  
  expect_equal(
    uneven_coverage_case$intron_start,
    11
  )
  
  expect_equal(
    uneven_coverage_case$intron_end,
    20
  ) 
  
})
