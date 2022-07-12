
# constructor
new_GenomicHistogram = function(x, interval_start = NULL, interval_end = NULL, region_id = NULL, chr = NULL, strand = NULL){

  # Checking types
  stopifnot(is.character(chr))
  stopifnot(is.character(strand))

  # Creating object
  new_Histogram(
    x,
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
  chr = attr(x, "chr")

  # Validate
  # 1. chr is of length 1
  if( length(chr) > 1 ) {
    stop("`chr` must have length 1", call. = FALSE)
  }

  # 2. validate_Histogram passes
  validate_Histogram(x)
}

# helper
#' Generates an S3 `GenomicHistogram` object
#'
#' @param x vector of counts/density
#' @param interval_start integer vector representing the starts of intervals
#' @param interval_end integer vector representing the ends of intervals
#' @param chr chromosome name
#' @param strand strand
#'
#' @return A GenomicHistogram object
#' @export
#'
#' @examples
#' x = GenomicHistogram(x = runif(10), interval_start = 1:10, interval_end = 1:10, chr = "chr1", strand = "+")
GenomicHistogram = function(x = double(), interval_start = integer(), interval_end = integer(), region_id = character(), chr = character(), strand = c("*", "+", "-")){

  # Coercing values to the right thing
  strand = match.arg(strand)

  if(length(x) > 0){
    if( missing(interval_start) & missing(interval_end)){
      interval_start = interval_end = seq(1, length(x), 1)
    } else if (missing(interval_start)){
      interval_start = interval_end
    } else if (missing(interval_end)){
      interval_end = interval_start
    }
  }

  if (!is.double(x)) {
    x = as.double(x)
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
  if( missing(region_id) & length(x) > 0 ){
    region_id = paste0(chr, ":", interval_start[1], "-", interval_end[length(x)], ":", strand)
  }
  if(!is.character(region_id)){
    region_id = as.character(region_id)
  }

  # Validate and return object
  validate_GenomicHistogram(new_GenomicHistogram(x = x, interval_start = interval_start, interval_end = interval_end, chr = chr, strand = strand, region_id = region_id))

}
