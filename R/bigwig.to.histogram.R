#' Imports BigWig file and generates histograms
#'
#' @param filename TODO
#' @param strand TODO
#' @param score.threshold TODO
#' @param regions TODO
#' @param gtf.file TODO
#' @param gene.or.transcript TODO
#' @param histogram.bin.size TODO
#' @param ... TODO
#'
#' @return
#' @export
bigwig.to.histogram = function(
  filename,
  strand  = c("+", "-", "."),
  score.threshold = 1,
  regions = NULL,
  gtf.file = NULL,
  gene.or.transcript = c("gene", "transcript"),
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
      gene.or.transcript = gene.or.transcript,
      select.strand = if(strand %in% c("+", "-")) strand else NULL,
      ...)
  } else {
    regions = GenomicRanges::reduce(bigwig.gr)
    ids = generate.identifiers(regions)
    regions = S4Vectors::split(regions, f = ids)
  }

  # Calculating Coverage
  histogram.coverage = coverage.to.histogram(
    coverage.rle = bigwig.coverage,
    regions = regions,
    histogram.bin.size = histogram.bin.size)

  # Returning list of histograms, region coordinates
  list(
    histogram.coverage = histogram.coverage,
    gene.model = regions,
    histogram.bin.size = histogram.bin.size)
}
