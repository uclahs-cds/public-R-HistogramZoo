context("Testing input and parameter variations to transcript.bed.to.histogram")

# Input Data
datadir = system.file("extdata", "rna_bedfiles",  package = "ConsensusPeaks")
filenames = file.path(datadir, paste0("Sample.", 1:20, ".bed"))
gtf.file = system.file("extdata", "genes.gtf",  package = "ConsensusPeaks")


test_that("Standard input yields the correct result format", {
  
  histograms = transcript.bed.to.histogram(
    filenames = filenames,
    n_fields = 12,
    gtf.file = gtf.file,
    gene.or.transcript = "gene",
    histogram.bin.size = 1
  )
  
  expect_length(histograms, 3)
  
  expect_named(
    histograms, 
    c("histogram.coverage","gene.model", "histogram.bin.size"))
  
  expect_length(histograms$gene.model, 2)
  expect_equal(
    names(histograms$gene.model),
    c( "ENSG00000178951.9", "ENSG00000185129.7")
  )
  
  expect_length(histograms$histogram.coverage, 2)
  expect_equal(
    names(histograms$histogram.coverage),
    c( "ENSG00000178951.9", "ENSG00000185129.7")    
  )
  
  expect_equal(
    histograms$histogram.bin.size,
    1
  )
  
})

test_that("Testing that selecting strand/transcript/chrs/ids return the correct results", {
  
  histograms.pos = transcript.bed.to.histogram(
    filenames = filenames,
    n_fields = 12,
    gtf.file = gtf.file,
    gene.or.transcript = "gene",
    histogram.bin.size = 1,
    select.strand = "+"
  )
  
  expect_named(
    histograms.pos$gene.model, 
    c("ENSG00000185129.7"))
  
  
  histograms.chr = transcript.bed.to.histogram(
    filenames = filenames,
    n_fields = 12,
    gtf.file = gtf.file,
    gene.or.transcript = "gene",
    histogram.bin.size = 1,
    select.chrs = "chr5"
  )
  
  expect_named(
    histograms.chr$gene.model,
    c("ENSG00000185129.7"))
  
  histograms.gene = transcript.bed.to.histogram(
    filenames = filenames,
    n_fields = 12,
    gtf.file = gtf.file,
    gene.or.transcript = "gene",
    histogram.bin.size = 1,
    select.ids = "ENSG00000178951.9"
  )
  
  expect_named(
    histograms.gene$gene.model,
    c("ENSG00000178951.9"))
  
  expect_warning(
    transcript.bed.to.histogram(
      filenames = filenames,
      n_fields = 12,
      gtf.file = gtf.file,
      gene.or.transcript = "gene",
      histogram.bin.size = 1,
      select.ids = "XX"
    )
  )
  
  expect_warning(
    transcript.bed.to.histogram(
      filenames = filenames,
      n_fields = 12,
      gtf.file = gtf.file,
      gene.or.transcript = "transcript",
      histogram.bin.size = 1
    )
  )
  
  
})

test_that("Testing that varying n_fields returns the correct results", {
  
  expect_error(
    transcript.bed.to.histogram(
      filenames = filenames,
      n_fields = 3,
      gtf.file = gtf.file,
      gene.or.transcript = "gene",
      histogram.bin.size = 1
    )
  )
  
  expect_error(
    transcript.bed.to.histogram(
      filenames = filenames,
      n_fields = 4,
      gtf.file = gtf.file,
      gene.or.transcript = "gene",
      histogram.bin.size = 1
    ),
    NA
  )
  
  expect_error(
    transcript.bed.to.histogram(
      filenames = filenames,
      n_fields = 6,
      gtf.file = gtf.file,
      gene.or.transcript = "gene",
      histogram.bin.size = 1
    ),
    NA
  )
  
})

test_that("Testing that varying bin size yields correct results", {
  
  histograms.5 = transcript.bed.to.histogram(
    filenames = filenames,
    n_fields = 12,
    gtf.file = gtf.file,
    gene.or.transcript = "gene",
    histogram.bin.size = 5
  )
  
  histograms.10 = transcript.bed.to.histogram(
    filenames = filenames,
    n_fields = 12,
    gtf.file = gtf.file,
    gene.or.transcript = "gene",
    histogram.bin.size = 10
  )

    # Creating bins
  bins.5 = GenomicRanges::tile(x = histograms.5$gene.model[[1]], width = 5)
  bins.5 = unlist(bins.5)
  bins.5 = GenomicRanges::width(bins.5)
  
  bins.10 = GenomicRanges::tile(x = histograms.10$gene.model[[1]], width = 10)
  bins.10 = unlist(bins.10)
  bins.10 = GenomicRanges::width(bins.10)
  
  # Checking that the binnedAverage add up to the same thing
  expect_equal(
    sum(bins.5*histograms.5$histogram.coverage[[1]]), 
    sum(bins.10*histograms.10$histogram.coverage[[1]]))
  
})