#' Calculates coverage of genes from annotated RNA bed files
#'
#' @param filenames A vector of BED filenames. The `name` column of the BED files must indicate gene or transcript name
#' @param n_fields Number of columns in the BED file that conform to BED file standards
#' @param gtf A GTF file
#' @param histogram_bin_size The bin size (base-pairs) to bin signal into a histogram
#' @param ... Additional parameters to be passed into GTF_to_GRangesList
#'
#' @return A list consisting of a list of histograms, a list of gene models and the histogram bin size
#'
#' @examples \dontrun{
#' file.directory = system.file("extdata", "rna_bedfiles", package = "ConsensusPeaks")
#' filenames = file.path(file.directory, paste0("Sample.", 1:20, ".bed"))
#' gtf = system.file("extdata", "genes.gtf", package = "ConsensusPeaks")
#'
#' histograms = transcript_BED_to_histogram(
#' filenames = filenames,
#' n_fields = 12,
#' gtf = gtf,
#' gene_or_transcript = "gene",
#' histogram_bin_size = 10)
#' }
#'
#' @export
transcript_BED_to_histogram = function(
  filenames,
  n_fields = c(4, 6, 12),
  gtf = NULL,
  histogram_bin_size = 1,
  gene_or_transcript = c("gene", "transcript"),
  select_strand = c("*", "+", "-"),
  select_chrs = NULL,
  select_ids = NULL,
  ...
){

  peaks = lapply(filenames, function(filename){
    segs = valr::read_bed( filename, n_fields = n_fields, comment = "#")
    if(n_fields == 12){ segs = valr::bed12_to_exons( segs ) }
    segs_gr = GenomicRanges::makeGRangesFromDataFrame( segs, keep.extra.columns = T )
    segs_gr = base0_to_base1(segs_gr)
    segs_gr = S4Vectors::split(segs_gr, f = segs_gr$name)
    segs_gr = GenomicRanges::reduce(segs_gr)
    unlist(segs_gr)
  })
  peaks = do.call(c, peaks)
  peaks = split(peaks, f = names(peaks))

  regions = GTF_to_GRangesList(
    gtf = gtf,
    gene_or_transcript = gene_or_transcript,
    select_strand = select_strand,
    select_chrs = select_chrs,
    select_ids = select_ids
  )

  ids = intersect(names(peaks), names(regions))
  if(length(ids) == 0) warning("No intersecting IDs between regions and peaks!")

  histogram_coverage =  vector("list", length(ids))
  names(histogram_coverage) = ids
  for(i in ids){
    peaks_cov = GenomicRanges::coverage(peaks[[i]])
    bins = GenomicRanges::tile(x = regions[[i]], width = histogram_bin_size)
    bins = unlist(bins)
    GenomeInfoDb::seqlevels(bins) = GenomeInfoDb::seqlevels(peaks_cov)
    cvg = GenomicRanges::binnedAverage(
      bins = bins,
      numvar = peaks_cov,
      varname = "cvg")
    histogram_coverage[[i]] <- new_GenomicHistogram(
      histogram_data = cvg$cvg,
      interval_start = GenomicRanges::start(cvg),
      interval_end = GenomicRanges::end(cvg),
      region_id = i,
      chr = as.character(GenomicRanges::seqnames(cvg))[1],
      strand = as.character(GenomicRanges::strand(cvg))[1]
    )
  }

  return(histogram_coverage)
}
