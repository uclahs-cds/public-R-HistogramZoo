context("Testing input and parameter variations to bigwig.to.histogram")

filename = system.file("extdata", "bigwigs", "test.bw", package = "ConsensusPeaks")
gtf = system.file("extdata", "genes.gtf",  package = "ConsensusPeaks")

strand = "."
score.threshold = 0
histogram.bin.size = 5


test_that("Pre-specified regions yield the correct result", {

  # Gene ENSG00000178951.9
  regions = GenomicRanges::GRanges(
    seqnames = "chr19",
    IRanges::IRanges(start = c(4043303),
                     end = c(4048244)),
    strand = "-")

  # With pre-specified regions
  histograms.gr.regions = bigwig.to.histogram(
    filename = filename,
    strand = strand,
    score.threshold = score.threshold,
    regions = regions,
    gtf.file = NULL,
    histogram.bin.size = histogram.bin.size)
  
  # Result format
  expect_length(histograms.gr.regions, 3)
  expect_named(
    histograms.gr.regions, 
    c("histogram.coverage","gene.model", "histogram.bin.size"))
  
  # Length of histograms
  expect_length(histograms.gr.regions$histogram.coverage, 1)
  
  # Number of bins
  expect_equal(
    length(histograms.gr.regions$histogram.coverage[[1]]), 
    ceiling(GenomicRanges::width(regions)[1]/histogram.bin.size))
  
  # Gene model and histograms length and names
  expect_equal(
    length(histograms.gr.regions$histogram.coverage), 
    length(histograms.gr.regions$gene.model))
  expect_equal(
    names(histograms.gr.regions$gene.model), 
    names(histograms.gr.regions$histogram.coverage))
  
  # Names of gene model
  expect_equal(
    names(histograms.gr.regions$gene.model)[1], 
    generate.identifiers(regions)[1])

  # Correct histogram bin size
  expect_equal(
    histogram.bin.size, 
    histograms.gr.regions$histogram.bin.size)
  
})

test_that("Testing that unspecified regions yield correct results", {
  
  # With no specified regions
  histograms.no.regions = bigwig.to.histogram(
    filename = filename,
    strand = strand,
    score.threshold = score.threshold,
    regions = NULL,
    gtf.file = NULL,
    histogram.bin.size = histogram.bin.size)

  # Result format
  expect_length(histograms.no.regions, 3)
  expect_named(
    histograms.no.regions, 
    c("histogram.coverage","gene.model", "histogram.bin.size"))
  
  # Length of histograms
  expect_length(histograms.no.regions$histogram.coverage, 2)
  
  # Selects only regions above 0
  expect_true(all(histograms.no.regions$histogram.coverage[[1]] > 0))
  
  # Gene model and histogram length and names
  expect_equal(
    length(histograms.no.regions$histogram.coverage), 
    length(histograms.no.regions$gene.model))
  expect_equal(
    names(histograms.no.regions$gene.model), 
    names(histograms.no.regions$histogram.coverage))
  
  # Names of gene model
  expect_named(
    histograms.no.regions$gene.model, 
    c("chr19:4043399-4043499","chr5:140114049-140114174")) 
  
  # Correct histogram bin size
  expect_equal(
    histogram.bin.size, 
    histograms.no.regions$histogram.bin.size)
  
})

test_that("Testing that importing a GTF file yields correct results", {

  # With a GTF file
  histograms.gtf = bigwig.to.histogram(
    filename = filename,
    strand = strand,
    score.threshold = score.threshold,
    regions = NULL,
    gtf.file = gtf,
    histogram.bin.size = histogram.bin.size)
  
  # Result format
  expect_length(histograms.gtf, 3)
  expect_named(
    histograms.gtf, 
    c("histogram.coverage","gene.model", "histogram.bin.size"))
  
  # Length of histograms
  expect_length(histograms.gtf$histogram.coverage, 2)
  
  # Gene model length and length of histograms
  expect_equal(
    length(histograms.gtf$histogram.coverage), 
    length(histograms.gtf$gene.model))
  expect_equal(
    names(histograms.gtf$gene.model), 
    names(histograms.gtf$histogram.coverage))
  
  # Takes the gene names
  expect_named(
    histograms.gtf$gene.model, 
    c("ENSG00000178951.9", "ENSG00000185129.7")) 
  
  # Correct histogram bin size
  expect_equal(
    histogram.bin.size, 
    histograms.gtf$histogram.bin.size)

})

test_that("Testing that varying bin size yields correct results", {
  
  regions = GenomicRanges::GRanges(
    seqnames = "chr19",
    IRanges::IRanges(start = c(4043303),
                     end = c(4048244)),
    strand = "-")
  
  histograms.5 = bigwig.to.histogram(
    filename = filename,
    strand = strand,
    score.threshold = score.threshold,
    regions = regions,
    gtf.file = NULL,
    histogram.bin.size = 5)
  
  expect_equal(
    length(histograms.5$histogram.coverage[[1]]), 
    ceiling(GenomicRanges::width(regions)[1]/5))

  histograms.10 = bigwig.to.histogram(
    filename = filename,
    strand = strand,
    score.threshold = score.threshold,
    regions = regions,
    gtf.file = NULL,
    histogram.bin.size = 10)
  
  expect_equal(
    length(histograms.10$histogram.coverage[[1]]), 
    ceiling(GenomicRanges::width(regions)[1]/10))
  
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

test_that("Testing that selecting for strand yields correct results", {

  histograms.pos = bigwig.to.histogram(
    filename = filename,
    strand = "+",
    score.threshold = score.threshold,
    regions = NULL,
    gtf.file = gtf,
    histogram.bin.size = histogram.bin.size)
  
  expect_named(
    histograms.pos$gene.model, 
    "ENSG00000185129.7") 
  
  histograms.neg = bigwig.to.histogram(
    filename = filename,
    strand = "-",
    score.threshold = score.threshold,
    regions = NULL,
    gtf.file = gtf,
    histogram.bin.size = histogram.bin.size)
  
  expect_named(
    histograms.neg$gene.model,
    "ENSG00000178951.9"
  )

  histograms.neutral = bigwig.to.histogram(
    filename = filename,
    strand = ".",
    score.threshold = score.threshold,
    regions = NULL,
    gtf.file = gtf,
    histogram.bin.size = histogram.bin.size)
  
  expect_named(
    histograms.neutral$gene.model,
    c("ENSG00000178951.9","ENSG00000185129.7")
  )

})

test_that("Testing that varying the score threshold yields correct results", {
  
  histograms.score.1 = bigwig.to.histogram(
    filename = filename,
    strand = strand,
    score.threshold = 1,
    regions = NULL,
    gtf.file = NULL,
    histogram.bin.size = histogram.bin.size)
  
  expect_true(all(histograms.score.1$histogram.coverage[[1]] > 1))
  
  histograms.score.2 = bigwig.to.histogram(
    filename = filename,
    strand = strand,
    score.threshold = 2,
    regions = NULL,
    gtf.file = NULL,
    histogram.bin.size = histogram.bin.size)
  
  expect_true(all(histograms.score.2$histogram.coverage[[1]] > 2))
  
})