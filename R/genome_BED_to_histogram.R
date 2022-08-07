
#' Generate Histograms from BED files
#'
#' @param filenames A vector of BED filenames. The `name` column of the BED files must indicate gene or transcript name
#' @param n_fields Number of columns in the BED file that conform to BED file standards
#' @param regions A GRanges or GRangesList object representing regions of interest defining histograms, a GRangesList object allows for the specification of non-continuous segments, default NULL - regions are taken as the set of elements with nonzero coverage
#' @param histogram_bin_size The bin size (base-pairs) to bin signal into a histogram
#' @param allow_overlapping_segments_per_sample logical, if FALSE, overlapping segments in the same file will be de-duplicated in the coverage calculation, if TRUE, they will be taken as separate input, default FALSE
#'
#' @return A list of GenomicHistogram objects
#' @export
genome_BED_to_histogram = function(
  filenames,
  n_fields = c(4, 6, 12),
  regions = NULL,
  histogram_bin_size = 1,
  allow_overlapping_segments_per_sample = F
){
  
  # Error checking
  existing_files <- file.exists(filenames)
  if(!all(existing_files)){
    stop(paste0(basename(filenames[!existing_files]), collapse = ","), " don't exist")
  }
  if(!(inherits(regions, "GRanges") | inherits(regions, "CompressedGRangesList") | is.null(regions))){
    stop("regions must be a GRanges or GRangesList object or set to NULL")
  }
  if(!is_equal_integer(histogram_bin_size) | histogram_bin_size < 1 ){
    stop("histogram_bin_size must be a positive integer")
  }
  if(!is.logical(allow_overlapping_segments_per_sample)){
    stop("allow_overlapping_segments_per_sample must be logical")
  }

  # Reading in peak files
  peaks <- lapply(filenames, function(filename){
    segs <- valr::read_bed( filename, n_fields = n_fields, comment = "#")
    if(n_fields == 12){ segs <- valr::bed12_to_exons( segs ) }
    segs_gr <- GenomicRanges::makeGRangesFromDataFrame( segs, keep.extra.columns = T )
    segs_gr <- base0_to_base1(segs_gr)
    if(!allow_overlapping_segments_per_sample){ segs_gr <- GenomicRanges::reduce(segs_gr) }
    return(segs_gr) 
  })
  peaks <- do.call(c, peaks)
  coverage <- GenomicRanges::coverage(peaks)

  # Formatting/generating regions, output: a GRangesList object defining regions
  if(inherits(regions, "CompressedGRangesList")){
    if(is.null(names(regions))){
      names(regions) <- generate_identifiers( unlist( range(regions) ) )
    }
  } else if (inherits(regions, "GRanges")){
    ids <- if(!is.null(regions$name)) regions$name else generate_identifiers(regions)
    regions <- S4Vectors::split(regions, f = ids)
  } else { # regions is NULL
    regions <- GenomicRanges::reduce(peaks)
    ids <- generate_identifiers(regions)
    regions <- S4Vectors::split(regions, f = ids)
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
