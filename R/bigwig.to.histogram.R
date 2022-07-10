#' Imports BigWig file and generates histograms
#'
#' @param filename Name of bigwig file for import
#' @param strand The strand of bigwig file from which the data originates. Default ".". If strand is "+" or "-", the strand will also be used to select regions of matching strand.
#' @param score.threshold A hard threshold for the score of the bigwig file. Scores higher than the threshold will be used in the computation of the histogram.
#' @param regions A GRanges object to select regions of interest
#' @param gtf.file A GTF file to select regions of interest
#' @param histogram.bin.size The bin size (base-pairs) to bin signal into a histogram
#' @param ... Additional parameters to be passed into gtf.to.genemodel
#'
#' @return A list consisting of a list of histograms, a list of gene models and the histogram bin size
#'
#' @examples \dontrun{
#' filename = system.file("extdata", "bigwigs",  "S1.bw", package = "ConsensusPeaks")
#' regions = GenomicRanges::GRanges(
#' seqnames = "chr1",
#' IRanges::IRanges(start = c(17950, 19350),
#'                  end = c(18000, 19600)),
#' strand = "*")
#'
#' histograms = BigWigToHistogram(
#' filename = filename,
#' regions = regions,
#' histogram.bin.size = 10)
#' }
#'
#'
#' @export
BigWigToHistogram = function(
  filename,
  strand  = c("*", "+", "-"),
  score.threshold = 0,
  regions = NULL,
  gtf.file = NULL,
  histogram.bin.size = 50,
  ...
){

  # Load BigWig file
  bigwig = valr::read_bigwig(filename, set_strand = strand)
  bigwig.gr = GenomicRanges::makeGRangesFromDataFrame(bigwig, keep.extra.columns = T)
  bigwig.gr = bigwig.gr[bigwig.gr$score > score.threshold]
  bigwig.coverage = GenomicRanges::coverage(bigwig.gr, weight = "score")

  # Loading regions
  if(!is.null(regions)){
    if(strand %in% c("+", "-")){
      regions = regions[GenomicRanges::strand(regions) == strand]
    }
    ids = if(!is.null(regions$name)) regions$name else generate.identifiers(regions)
    regions = S4Vectors::split(regions, f = ids)
  } else if(!is.null(gtf.file)){
    regions = gtf.to.genemodel(
      gtf.file = gtf.file,
      select.strand = if(strand %in% c("+", "-")) strand else NULL,
      ...)
  } else {
    regions = GenomicRanges::reduce(bigwig.gr)
    ids = generate.identifiers(regions)
    regions = S4Vectors::split(regions, f = ids)
  }

  # Calculating Coverage
  coverage.to.histogram(
    coverage.rle = bigwig.coverage,
    regions = regions,
    histogram.bin.size = histogram.bin.size
  )
}
