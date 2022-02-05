#' Calculates coverage of genes from annotated RNA bed files
#'
#' @param filenames TODO
#' @param n_fields TODO
#' @param regions.of.interest TODO
#' @param gtf.file TODO
#' @param gene.or.transcript TODO
#' @param histogram.bin.size TODO
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
transcript.bed.to.hist = function(
  filenames,
  n_fields = 3,
  regions.of.interest = NULL,
  gtf.file = NULL,
  gene.or.transcript = c("gene", "transcript"),
  histogram.bin.size = 1,
  ...
){

  peaks = lapply(filenames, function(filename){
    segs = valr::read_bed( filename, n_fields = n_fields, ...)
    if(n_fields == 12){ segs = valr::bed12_to_exons( segs ) }
    segs.gr = GenomicRanges::makeGRangesFromDataFrame( segs, keep.extra.columns = T )
    segs.gr = base0.to.base1(segs.gr)
    segs.gr = S4Vectors::split(segs.gr, f = segs.gr$name)
    segs.gr = GenomicRanges::reduce(segs.gr)
    unlist(segs.gr)
  })
  peaks = do.call(c, peaks)
  peaks = split(peaks, f = names(peaks))

  genes = names(peaks)
  if(!is.null(regions.of.interest)){ genes = intersect(genes, regions.of.interest)}

  regions = gtf.to.genemodel(
    gtf.file = gtf.file,
    split.by.strand = F,
    gene.or.transcript = gene.or.transcript,
    select.chrs = NULL,
    select.genes = genes)

  # TODO: Error check to make sure names and genes match

  histogram.coverage =  vector("list", length(genes))
  names(histogram.coverage) = genes
  for(i in genes){
    peaks.cov = GenomicRanges::coverage(peaks[[i]])
    bins = GenomicRanges::tile(x = regions[[i]], width = histogram.bin.size)
    bins = unlist(bins)
    GenomeInfoDb::seqlevels(bins) = GenomeInfoDb::seqlevels(peaks.cov)
    cvg = GenomicRanges::binnedAverage(
      bins = bins,
      numvar = peaks.cov,
      varname = "cvg")
    cvg = cvg$cvg
    names(cvg) = NULL
    histogram.coverage[[i]] <- cvg
  }

  list(
    histogram.coverage = histogram.coverage,
    gene.model = regions,
    histogram.bin.size = histogram.bin.size)
}
