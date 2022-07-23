
# constructor
new_GenomicHistogram = function(histogram_data = NULL, interval_start = NULL, interval_end = NULL, region_id = NULL, chr = NULL, strand = NULL){

  # Checking types
  stopifnot(is.character(chr))
  stopifnot(is.character(strand))

  # Creating object
  new_Histogram(
    histogram_data = histogram_data,
    interval_start = interval_start,
    interval_end = interval_end,
    region_id = region_id,
    chr = chr,
    strand = strand,
    class = "GenomicHistogram"
  )
}

# validator
validate_GenomicHistogram = function(x){

  # Attributes
  chr = x$chr

  # Validate
  # 1. chr is of length 1
  if( length(chr) > 1 ) {
    stop("chr must have length 1", call. = FALSE)
  }

  # 2. validate_Histogram passes
  validate_Histogram(x)
}

# helper
#' Generates an S3 `GenomicHistogram` object
#'
#' @param histogram_data vector of counts/density
#' @param interval_start integer vector representing the starts of intervals
#' @param interval_end integer vector representing the ends of intervals
#' @param chr chromosome name
#' @param strand strand
#'
#' @return A GenomicHistogram object
#' @export
#'
#' @examples
#' x = GenomicHistogram(histogram_data = runif(10), interval_start = 1:10, interval_end = 1:10, chr = "chr1", strand = "+")
GenomicHistogram = function(histogram_data = double(), interval_start = integer(), interval_end = integer(), region_id = character(), chr = character(), strand = c("*", "+", "-")){

  # Coercing values to the right thing
  strand = match.arg(strand)

  if(length(histogram_data) > 0){
    if( missing(interval_start) & missing(interval_end)){
      interval_start = interval_end = seq(1, length(histogram_data), 1)
    } else if (missing(interval_start)){
      interval_start = interval_end
    } else if (missing(interval_end)){
      interval_end = interval_start
    }
  }

  if (!is.double(histogram_data)) {
    histogram_data = as.double(histogram_data)
  }
  if (!is.integer(interval_start)){
    interval_start = as.integer(interval_start)
  }
  if (!is.integer(interval_end)){
    interval_end = as.integer(interval_end)
  }
  if(!is.character(chr)){
    chr = as.character(chr)
  }

  # Region ID
  if( missing(region_id) & length(histogram_data) > 0 ){
    region_id = paste0(chr, ":", interval_start[1], "-", interval_end[length(histogram_data)], ":", strand)
  }
  if(!is.character(region_id)){
    region_id = as.character(region_id)
  }

  # Validate and return object
  validate_GenomicHistogram(new_GenomicHistogram(histogram_data = histogram_data, interval_start = interval_start, interval_end = interval_end, chr = chr, strand = strand, region_id = region_id))

}

#' Title
#'
#' @param x
#' @param i
#'
#' @return
#' @export
#'
#' @examples
`[.GenomicHistogram` = function(x, i){
  new_GenomicHistogram(
    histogram_data = x$histogram_data[i],
    interval_start = x$interval_start[i],
    interval_end = x$interval_end[i],
    region_id = x$region_id,
    chr = x$chr,
    strand = x$strand)
}

reassign_region_id.GenomicHistogram = function(x, region_id){

  stopifnot(inherits(x, "Histogram"))

  # Creating a region id
  if(missing(region_id) & length(x$histogram_data) > 0){
    region_id = paste0(x$chr, ":", x$interval_start[1], "-", x$interval_end[length(x$histogramm_data)], ":", x$strand)
  }
  if(!is.character(region_id)){
    region_id = as.character(region_id)
  }

  x$region_id <- region_id

  return(x)
}
