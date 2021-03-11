#' Plots the gene coverage for a given gene
#'
#' @param PEAKS A data frame containing the following columns, and potentially extras, usually found in a BED12 file, base 0 system
#' @param GENES A character vector of genes to be tested
#' @param ANNOTATION Annotations returned from \link[ConsensusPeaks]{read.gtf}.
#'
#' @return ggplot of the gene coverage histogram
#' @export
#'
#' @examples
#' gtf.path = system.file("extdata/test.gtf", package = "ConsensusPeaks")
#' PARAMETERS = list("GTF" = gtf, "GENE" = "ENSGXX")
#' annotation = read.gtf(PARAMETERS)
#' ensgxx = simulate.gaussian.peaks(
#'   MU = c(100, 150),
#'   SD = c(10, 20),
#'   EXTEND.WIDTH = c(50, 25),
#'   NSAMPLES = c(10, 30),
#'   GENE = "ENSGXX",
#'   GTF = gtf,
#'   ANNOTATION = annotation,
#'   SEED = 123)
#' gene.coverage.hist(ensgxx, "ENSGXX", annotation)
gene.coverage.hist <- function(PEAKS, GENE, ANNOTATION) {

  # PEAKSGR
  PEAKSGR = .retrieve.peaks.as.granges(PEAKS = PEAKS, GENE = GENE, DF = F)

  # Get Gene Information
  GENEINFO = .get.gene.anno(GENE, ANNOTATION)

  # Converting to RNA
  GENEPEAKSGR = GenomicRanges::shift(PEAKSGR, -1*GENEINFO$left+1)
  GenomicRanges::start(GENEPEAKSGR) = GENEINFO$DNA2RNA[GenomicRanges::start(GENEPEAKSGR)+1]
  GenomicRanges::end(GENEPEAKSGR) = GENEINFO$DNA2RNA[GenomicRanges::end(GENEPEAKSGR)]

  # Reduce Overlapping Peaks in the Same Sample
  GENEPEAKSGR = S4Vectors::split(GENEPEAKSGR, GENEPEAKSGR$sample)
  GENEPEAKSGR = unlist(GenomicRanges::reduce(GENEPEAKSGR))

  REDUCED.GENE.PEAKS.GR = GenomicRanges::reduce(GENEPEAKSGR)
  PEAK.COVERAGE = GenomicRanges::coverage(GENEPEAKSGR)
  BINS = unlist(GenomicRanges::tile(REDUCED.GENE.PEAKS.GR, width = 1))
  BIN.COUNTS = data.frame(GenomicRanges::binnedAverage(BINS, PEAK.COVERAGE, "Coverage"), stringsAsFactors = F)

  ggplot2::ggplot() +
    ggplot2::geom_step(data=BIN.COUNTS, ggplot2::aes(x=start, y=Coverage), colour='grey') +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(GENE) +
    ggplot2::ylab("Coverage (at BP resolution)") +
    ggplot2::xlab("Transcript Coordinate")
}
