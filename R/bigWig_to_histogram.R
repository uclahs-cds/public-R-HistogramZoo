#' Imports a single bigWig file and generates histograms
#'
#' @param filename Name of bigwig file for import
#' @param strand The strand of bigwig file from which the data originates. Default ".". If strand is "+" or "-", the strand will also be used to select regions of matching strand.
#' @param score_threshold A hard threshold for the score of the bigwig file. Scores higher than the threshold will be used in the computation of the histogram.
#' @param regions A GRanges object to select regions of interest
#' @param gtf A GTF file to select regions of interest
#' @param histogram_bin_size The bin size (base-pairs) to bin signal into a histogram
#' @param ... Additional parameters to be passed into GTF_to_GRangesList
#'
#' @return A list consisting of GenomicHistograms
#'
#' @examples \dontrun{
#' filename = system.file("extdata", "bigwigs",  "S1.bw", package = "ConsensusPeaks")
#' regions = GenomicRanges::GRanges(
#' seqnames = "chr1",
#' IRanges::IRanges(start = c(17950, 19350),
#'                  end = c(18000, 19600)),
#' strand = "*")
#'
#' histograms = bigWig_to_histogram(
#' filename = filename,
#' regions = regions,
#' histogram_bin_size = 10)
#' }
#'
#'
#' @export
bigWig_to_histogram = function(
  filename,
  strand  = c("*", "+", "-"),
  score_threshold = 0,
  regions = NULL,
  gtf = NULL,
  histogram_bin_size = 50,
  ...
){

  # Load BigWig file
  bigwig = valr::read_bigwig(filename, set_strand = strand)
  bigwig_gr = GenomicRanges::makeGRangesFromDataFrame(bigwig, keep.extra.columns = T)
  bigwig_gr = bigwig_gr[bigwig_gr$score > score_threshold]
  bigwig_coverage = GenomicRanges::coverage(bigwig_gr, weight = "score")

  # Loading regions
  if(!is.null(regions)){
    if(strand %in% c("+", "-")){
      regions = regions[GenomicRanges::strand(regions) == strand]
    }
    ids = if(!is.null(regions$name)) regions$name else generate_identifiers(regions)
    regions = S4Vectors::split(regions, f = ids)
  } else if(!is.null(gtf)){
    regions = GTF_to_GRangesList(
      gtf = gtf,
      select_strand = if(strand %in% c("+", "-")) strand else NULL,
      ...)
  } else {
    regions = GenomicRanges::reduce(bigwig_gr)
    ids = generate_identifiers(regions)
    regions = S4Vectors::split(regions, f = ids)
  }

  # Calculating Coverage
  return(
    coverage_to_histogram(
      coverage = bigwig_coverage,
      regions = regions,
      histogram_bin_size = histogram_bin_size
    )
  )
}
