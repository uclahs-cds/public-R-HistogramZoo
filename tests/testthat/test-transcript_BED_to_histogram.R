context("Testing input and parameter variations to transcript_BED_to_histogram")

# Input Data
datadir = system.file("extdata", "rna_bedfiles",  package = "ConsensusPeaks")
filenames = file.path(datadir, paste0("Sample.", 1:20, ".bed"))
gtf = system.file("extdata", "genes.gtf",  package = "ConsensusPeaks")


test_that("Standard input yields the correct result format", {

  histograms = transcript_BED_to_histogram(
    filenames = filenames,
    n_fields = 12,
    gtf = gtf,
    gene_or_transcript = "gene",
    histogram_bin_size = 1
  )

  expect_length(histograms, 2)

  expect_named(
    histograms,
    c("ENSG00000178951.9", "ENSG00000185129.7"))

  # TODO: test that each is a valid histogram object

})

test_that("Testing that selecting strand/transcript/chrs/ids return the correct results", {

  histograms.pos = transcript_BED_to_histogram(
    filenames = filenames,
    n_fields = 12,
    gtf = gtf,
    gene_or_transcript = "gene",
    histogram_bin_size = 1,
    select_strand = "+"
  )

  expect_named(
    histograms.pos,
    "ENSG00000185129.7")


  histograms.chr = transcript_BED_to_histogram(
    filenames = filenames,
    n_fields = 12,
    gtf = gtf,
    gene_or_transcript = "gene",
    histogram_bin_size = 1,
    select_chrs = "chr5"
  )

  expect_named(
    histograms.chr,
    "ENSG00000185129.7")

  histograms.gene = transcript_BED_to_histogram(
    filenames = filenames,
    n_fields = 12,
    gtf = gtf,
    gene_or_transcript = "gene",
    histogram_bin_size = 1,
    select_ids = "ENSG00000178951.9"
  )

  expect_named(
    histograms.gene,
    "ENSG00000178951.9")

  expect_warning(
    transcript_BED_to_histogram(
      filenames = filenames,
      n_fields = 12,
      gtf = gtf,
      gene_or_transcript = "gene",
      histogram_bin_size = 1,
      select_ids = "XX"
    )
  )

  expect_warning(
    transcript_BED_to_histogram(
      filenames = filenames,
      n_fields = 12,
      gtf = gtf,
      gene_or_transcript = "transcript",
      histogram_bin_size = 1
    )
  )

})

test_that("Testing that varying n_fields returns the correct results", {

  expect_error(
    transcript_BED_to_histogram(
      filenames = filenames,
      n_fields = 3,
      gtf = gtf,
      gene_or_transcript = "gene",
      histogram_bin_size = 1
    )
  )

  expect_error(
    transcript_BED_to_histogram(
      filenames = filenames,
      n_fields = 4,
      gtf = gtf,
      gene_or_transcript = "gene",
      histogram_bin_size = 1
    ),
    NA
  )

  expect_error(
    transcript_BED_to_histogram(
      filenames = filenames,
      n_fields = 6,
      gtf = gtf,
      gene_or_transcript = "gene",
      histogram_bin_size = 1
    ),
    NA
  )

})

test_that("Testing that varying bin size yields correct results", {

  histograms.5 = transcript_BED_to_histogram(
    filenames = filenames,
    n_fields = 12,
    gtf = gtf,
    gene_or_transcript = "gene",
    histogram_bin_size = 5
  )

  histograms.10 = transcript_BED_to_histogram(
    filenames = filenames,
    n_fields = 12,
    gtf = gtf,
    gene_or_transcript = "gene",
    histogram_bin_size = 10
  )

  # Creating bins
  bins.5 = histograms.5[[1]]$interval_end - histograms.5[[1]]$interval_start + 1
  bins.10 = histograms.10[[1]]$interval_end - histograms.10[[1]]$interval_start + 1

  # Checking that the binnedAverage add up to the same thing
  expect_equal(
    sum(bins.5*histograms.5[[1]]$histogram_data),
    sum(bins.10*histograms.10[[1]]$histogram_data)
  )

})
