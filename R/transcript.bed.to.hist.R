#' Calculates coverage of genes from annotated RNA bed files
#'
#' @param filenames A vector of BED filenames
#' @param n_fields Number of columns in the BED file that conform to BED file standards
#' @param gtf.file A GTF file
#' @param gene.or.transcript Whether histograms should be computed on gene annotations or transcript annotations
#' @param select.ids Select elements by matching ids to genes or transcripts (depending on gene.or.transcript)
#' @param histogram.bin.size The bin size (base-pairs) to bin signal into a histogram
#' @param ... Additional parameters to be passed into gtf.to.genemodel
#'
#' @return A list consisting of a list of histograms, a list of gene models and the histogram bin size
#'
#' @example \dontrun{
#' file.directory = system.file("extdata", "RNA_bedfiles", package = "ConsensusPeaks")
#' filenames = list.files(file.directory)
#' gtf.file = system.file("extdata", "genes.gtf", package = "ConsensusPeaks")
#'
#' histograms = transcript.bed.to.histogram(
#' filenames = filenames,
#' n_fields 12,
#' gtf.file = gtf.file,
#' gene.or.transcript = "gene",
#' histogram.bin.size = 10)
#' }
#'
#' @export
transcript.bed.to.histogram = function(
  filenames,
  n_fields = c(3, 4, 6, 12),
  gtf.file = NULL,
  gene.or.transcript = c("gene", "transcript"),
  select.ids = NULL,
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
