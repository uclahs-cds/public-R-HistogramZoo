context("genome_BED_to_histogram")

library(GenomicRanges)

datadir = system.file("extdata", "dna_bedfiles",  package = "ConsensusPeaks")
filenames = list.files(datadir, pattern = ".bed$")
filenames = file.path(datadir, filenames)
n_fields = 6

# Defining regions as a GRanges object
gr_regions = GRanges(
  seqnames = "chr9",
  IRanges(
    start = c(70417946, 70420005),
    end = c(70420000, 70421544)
  ),
  strand = "*"
)

# Defining regions as a GRangesList object
glist_regions = split(gr_regions, f = 1:2)

# GRangesList representes a set of segmented regions
glist_segmented = GRangesList(gr_regions)

# Test overlap
datadir_overlap = system.file("extdata", "dna_bed_test_overlap",  package = "ConsensusPeaks")
filenames_overlap = list.files(datadir_overlap, pattern = ".bed$")
filenames_overlap = file.path(datadir_overlap, filenames_overlap)

# Test strand
datadir_strand = system.file("extdata", "dna_bed_test_strand",  package = "ConsensusPeaks")
filenames_strand = list.files(datadir_strand, pattern = ".bed$")
filenames_strand = file.path(datadir_strand, filenames_strand)

test_that("Basic use", {
  
  histograms = genome_BED_to_histogram(
    filenames = filenames,
    n_fields = n_fields,
    histogram_bin_size = 1
  )
  
  expect_length(histograms[[1]], 3599)
  
  expect_equal(
    histograms[[1]]$interval_start, 
    histograms[[1]]$interval_end
  )
  
  expect_equal(
    histograms[[1]]$chr,
    "chr9"
  )
  
  expect_equal(
    histograms[[1]]$strand,
    "*"
  )
  
  expect_equal(
    histograms[[1]]$region_id,
    "chr9:70417946-70421544:*"
  )

})

test_that("User-defined regions", {
  
  # GRanges regions
  histograms_gr = genome_BED_to_histogram(
    filenames = filenames,
    n_fields = n_fields,
    histogram_bin_size = 1,
    regions = gr_regions
  )
  
  # GRangesList regions
  histograms_glist = genome_BED_to_histogram(
    filenames = filenames,
    n_fields = n_fields,
    histogram_bin_size = 1,
    regions = glist_regions
  )
  
  expect_length(histograms_gr, 2)
  expect_length(histograms_glist, 2)
  
  expect_named(
    histograms_gr,
    c("chr9:70417946-70420000:*", "chr9:70420005-70421544:*"))
  expect_named(
    histograms_glist, 
    c("1", "2"))
  
  expect_equal(
    histograms_gr[[1]]$interval_start,
    histograms_glist[[1]]$interval_start
  )
  
  expect_equal(
    histograms_gr[[1]]$histogram_data,
    histograms_glist[[1]]$histogram_data
  )
  
})

test_that("GRangesList works on segmented regions", {
  
  # GRangesList regions
  histograms_glist_segmented = genome_BED_to_histogram(
    filenames = filenames,
    n_fields = n_fields,
    histogram_bin_size = 1,
    regions = glist_segmented
  )
  
  expect_length(histograms_glist_segmented, 1)
  expect_named(histograms_glist_segmented, "chr9:70417946-70421544:*")
  expect_equal(
    setdiff(70417946:70421544, 70420001:70420004),
    histograms_glist_segmented[[1]]$interval_start)
  
})


test_that("Histogram bin size", {
  
  histograms = genome_BED_to_histogram(
    filenames = filenames,
    n_fields = n_fields,
    histogram_bin_size = 1
  )
  
  histograms_10 = genome_BED_to_histogram(
    filenames = filenames,
    n_fields = n_fields,
    histogram_bin_size = 10
  )
  
  expect_length(
    histograms_10[[1]], 
    ceiling(length(histograms[[1]])/10))
  
  expect_equal(
    sum(histograms[[1]]$histogram_data),
    sum(histograms_10[[1]]$histogram_data*(histograms_10[[1]]$interval_end - histograms_10[[1]]$interval_start + 1))
  )
})

test_that("Non-overlapping regions per file", {
  
  histogram_dup = genome_BED_to_histogram(
    filenames = filenames_overlap,
    n_fields = n_fields,
    histogram_bin_size = 1,
    allow_overlapping_segments_per_sample = T
  )
  
  histogram_unique = genome_BED_to_histogram(
    filenames = filenames_overlap,
    n_fields = n_fields,
    histogram_bin_size = 1,
    allow_overlapping_segments_per_sample = F
  )
  
  expect_equal(
    histogram_unique[[1]]$histogram_data*2,
    histogram_dup[[1]]$histogram_data
  )
  
})

test_that("Potential strand-based issues", {
  
  histograms_strand = genome_BED_to_histogram(
    filenames = filenames_strand,
    n_fields = n_fields,
    histogram_bin_size = 1
  )
  
  expect_length(histograms_strand, 3)
  expect_named(
    histograms_strand,
    c("chr9:70418161-70418463:*", "chr9:70418161-70418463:+", "chr9:70418161-70418463:-"))
  
})