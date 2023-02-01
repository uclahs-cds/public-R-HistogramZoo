context("bigWig_to_histogram")

filename <- system.file("extdata", "bigwigs", "test.bw", package = "HistogramZoo")
gtf <- system.file("extdata", "genes.gtf",  package = "HistogramZoo")

strand <- "*"
score_threshold <- 0
histogram_bin_size <- 5


test_that("Pre-specified regions yield the correct result", {

  # Gene ENSG00000178951.9
  regions <- GenomicRanges::GRanges(
    seqnames = "chr19",
    IRanges::IRanges(start = c(4043303),
                     end = c(4048244)),
    strand = "-")

  # With pre-specified regions
  histograms_gr_regions <- bigWig_to_histogram(
    filename = filename,
    strand = strand,
    score_threshold = score_threshold,
    regions = regions,
    gtf = NULL,
    histogram_bin_size = histogram_bin_size)

  # Result format
  expect_length(histograms_gr_regions, 1)
  expect_named(
    histograms_gr_regions,
    "chr19:4043303-4048244:-")

  # Number of bins
  expect_equal(
    length(histograms_gr_regions[[1]]$histogram_data),
    ceiling(GenomicRanges::width(regions)[1]/histogram_bin_size))

  # TODO: Check other aspects of the return object are correct

})

test_that("Testing that unspecified regions yield correct results", {

  # With no specified regions
  histograms_no_regions <- bigWig_to_histogram(
    filename = filename,
    strand = strand,
    score_threshold = score_threshold,
    regions = NULL,
    gtf = NULL,
    histogram_bin_size = histogram_bin_size)

  # Result format
  expect_length(histograms_no_regions, 2)
  expect_named(
    histograms_no_regions,
    c("chr19:4043399-4043499:*", "chr5:140114049-140114174:*"))

  # Selects only regions above 0
  expect_true(all(histograms_no_regions[[1]]$histogram_data > 0))

})

test_that("Testing that importing a GTF file yields correct results", {

  # With a GTF file
  histograms_gtf <- bigWig_to_histogram(
    filename = filename,
    strand = strand,
    score_threshold = score_threshold,
    regions = NULL,
    gtf = gtf,
    histogram_bin_size = histogram_bin_size)

  # Result format
  expect_length(histograms_gtf, 2)
  expect_named(
    histograms_gtf,
    c("ENSG00000178951.9","ENSG00000185129.7"))

})

test_that("Testing that varying bin size yields correct results", {

  regions <- GenomicRanges::GRanges(
    seqnames = "chr19",
    IRanges::IRanges(start = c(4043303),
                     end = c(4048244)),
    strand = "-")

  histograms_5 <- bigWig_to_histogram(
    filename = filename,
    strand = strand,
    score_threshold = score_threshold,
    regions = regions,
    gtf = NULL,
    histogram_bin_size = 5)

  expect_equal(
    length(histograms_5[[1]]$histogram_data),
    ceiling(GenomicRanges::width(regions)[1]/5))

  histograms_10 <- bigWig_to_histogram(
    filename = filename,
    strand = strand,
    score_threshold = score_threshold,
    regions = regions,
    gtf = NULL,
    histogram_bin_size = 10)

  expect_equal(
    length(histograms_10[[1]]$histogram_data),
    ceiling(GenomicRanges::width(regions)[1]/10))

  # Creating bins
  bins_5 <- histograms_5[[1]]$interval_end - histograms_5[[1]]$interval_start + 1
  bins_10 <- histograms_10[[1]]$interval_end - histograms_10[[1]]$interval_start + 1

  # Checking that the binnedAverage add up to the same thing
  expect_equal(
    sum(bins_5*histograms_5[[1]]$histogram_data),
    sum(bins_10*histograms_10[[1]]$histogram_data))

})

test_that("Testing that selecting for strand yields correct results", {

  histograms_pos <- bigWig_to_histogram(
    filename = filename,
    strand = "+",
    score_threshold = score_threshold,
    regions = NULL,
    gtf = gtf,
    histogram_bin_size = histogram_bin_size)

  expect_named(
    histograms_pos,
    "ENSG00000185129.7")

  histograms_neg <- bigWig_to_histogram(
    filename = filename,
    strand = "-",
    score_threshold = score_threshold,
    regions = NULL,
    gtf = gtf,
    histogram_bin_size = histogram_bin_size)

  expect_named(
    histograms_neg,
    "ENSG00000178951.9"
  )

  histograms_neutral <- bigWig_to_histogram(
    filename = filename,
    strand = "*",
    score_threshold = score_threshold,
    regions = NULL,
    gtf = gtf,
    histogram_bin_size = histogram_bin_size)

  expect_named(
    histograms_neutral,
    c("ENSG00000178951.9","ENSG00000185129.7")
  )

})

test_that("Testing that varying the score threshold yields correct results", {

  histograms_score_1 <- bigWig_to_histogram(
    filename = filename,
    strand = strand,
    score_threshold = 1,
    regions = NULL,
    gtf = NULL,
    histogram_bin_size = histogram_bin_size)

  expect_true(all(histograms_score_1[[1]]$histogram_data > 1))

  histograms_score_2 <- bigWig_to_histogram(
    filename = filename,
    strand = strand,
    score_threshold = 2,
    regions = NULL,
    gtf = NULL,
    histogram_bin_size = histogram_bin_size)

  expect_true(all(histograms_score_2[[1]]$histogram_data > 2))

})
