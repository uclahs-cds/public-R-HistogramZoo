
#' Title
#'
#' @param filenames
#' @param n_fields
#' @param regions
#' @param histogram_bin_size
#' @param gtf
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
genome_BED_to_histogram = function(
  filenames,
  n_fields = c(4, 6, 12),
  regions = NULL,
  histogram_bin_size = 1,
  gtf = NULL,
  ...
){

  peaks = lapply(filenames, function(filename){
    segs = valr::read_bed( filename, n_fields = n_fields, comment = "#")
    if(n_fields == 12){ segs = valr::bed12_to_exons( segs ) }
    segs.gr = GenomicRanges::makeGRangesFromDataFrame( segs, keep.extra.columns = T )
    segs.gr = base0_to_base1(segs.gr)
    segs.gr = GenomicRanges::reduce(segs.gr)
    segs.gr
  })
  peaks.gr = do.call(c, peaks)
  coverage = GenomicRanges::coverage(peaks.gr)

  # Loading regions
  if(!is.null(regions)){
    ids = if(!is.null(regions$name)) regions$name else generate_identifiers(regions)
    regions = S4Vectors::split(regions, f = ids)
  } else if(!is.null(gtf)){
    regions = GTF_to_GRangesList(gtf = gtf, ...)
  } else {
    regions = GenomicRanges::reduce(peaks.gr)
    ids = generate_identifiers(regions)
    regions = S4Vectors::split(regions, f = ids)
  }

  # Calculating Coverage
  return(
    coverage_to_histogram(
      coverage = coverage,
      regions = regions,
      histogram_bin_size = histogram_bin_size
    )
  )
}
