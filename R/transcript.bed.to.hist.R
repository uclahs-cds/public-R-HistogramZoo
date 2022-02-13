#' Calculates coverage of genes from annotated RNA bed files
#'
#' @param filenames TODO
#' @param n_fields TODO
#' @param gtf.file TODO
#' @param gene.or.transcript TODO
#' @param select.ids TODO
#' @param histogram.bin.size TODO
#' @param ... TODO
#'
#' @return
#' @export
transcript.bed.to.histogram = function(
  filenames,
  n_fields = 3,
  gtf.file = NULL,
  gene.or.transcript = c("gene", "transcript"),
  select.ids = NULL,
  histogram.bin.size = 1,
  ...
){

  peaks = lapply(filenames, function(filename){
    segs = valr::read_bed( filename, n_fields = n_fields, comment = "#")
    if(n_fields == 12){ segs = valr::bed12_to_exons( segs ) }
    segs.gr = GenomicRanges::makeGRangesFromDataFrame( segs, keep.extra.columns = T )
    segs.gr = ConsensusPeaks:::base0.to.base1(segs.gr)
    segs.gr = S4Vectors::split(segs.gr, f = segs.gr$name)
    segs.gr = GenomicRanges::reduce(segs.gr)
    unlist(segs.gr)
  })
  peaks = do.call(c, peaks)
  peaks = split(peaks, f = names(peaks))

  ids = names(peaks)
  if(!is.null(select.ids)){
    ids = intersect(ids, select.ids)
  }

  regions = gtf.to.genemodel(
    gtf.file = gtf.file,
    gene.or.transcript = gene.or.transcript,
    select.ids = ids,
    ...)

  histogram.coverage =  vector("list", length(ids))
  names(histogram.coverage) = ids
  for(i in ids){
    peaks.cov = GenomicRanges::coverage(peaks[[i]])
    bins = GenomicRanges::tile(x = regions[[i]], width = histogram.bin.size)
    bins = unlist(bins)
    GenomeInfoDb::seqlevels(bins) = GenomeInfoDb::seqlevels(peaks.cov)
    cvg = GenomicRanges::binnedAverage(
      bins = bins,
      numvar = peaks.cov,
      varname = "cvg")
    histogram.coverage[[i]] <- cvg$cvg
  }

  list(
    histogram.coverage = histogram.coverage,
    gene.model = regions,
    histogram.bin.size = histogram.bin.size)
}
