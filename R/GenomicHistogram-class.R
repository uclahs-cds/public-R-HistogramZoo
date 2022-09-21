
# constructor
#' GenomicHistogram constructor
#'
#' @param histogram_data vector of counts/density
#' @param interval_start integer vector representing the starts of intervals
#' @param interval_end integer vector representing the ends of intervals
#' @param region_id character identifier for the region of interest
#' @param chr chromosome name
#' @param strand strand
#'
#' @return A GenomicHistogram object
new_GenomicHistogram <- function(
    histogram_data = NULL,
    interval_start = NULL,
    interval_end = NULL,
    region_id = NULL,
    bin_width = 1L,
    chr = NULL,
    strand = NULL
  ){

  # Checking types
  stopifnot(is.character(chr))
  stopifnot(is.character(strand))

  # Creating object
  return(
    new_Histogram(
      histogram_data = histogram_data,
      interval_start = interval_start,
      interval_end = interval_end,
      region_id = region_id,
      bin_width = bin_width,
      chr = chr,
      strand = strand,
      class = "GenomicHistogram"
    )
  )
}

# validator
#' Validates GenomicHistogram objects
#'
#' @param x a GenomicHistogram object candidate
#' @return a validated GenomicHistogram object
validate_GenomicHistogram <- function(x){

  # Attributes
  chr <- x$chr

  # Validate
  # 1. chr is of length 1
  if( length(chr) > 1 ) {
    stop("chr must have length 1", call. = FALSE)
  }

  # 2. validate_Histogram passes
  return(
    validate_Histogram(x)
  )
}

# helper
#' Generates an S3 `GenomicHistogram` object
#'
#' @param histogram_data vector of counts/density
#' @param interval_start integer vector representing the starts of intervals
#' @param interval_end integer vector representing the ends of intervals
#' @param region_id character identifier for the region of interest
#' @param chr chromosome name
#' @param strand strand
#'
#' @return a GenomicHistogram object
#' @export
#'
#' @examples
#' x = GenomicHistogram(
#' histogram_data = runif(10),
#' interval_start = 1:10,
#' interval_end = 1:10,
#' chr = "chr1",
#' strand = "+")
GenomicHistogram <- function(
    histogram_data = double(),
    interval_start = integer(),
    interval_end = integer(),
    region_id = character(),
    chr = character(),
    strand = c("*", "+", "-")
){

  # Coercing values to the right thing
  strand <- match.arg(strand)

  if(length(histogram_data) > 0){
    if( missing(interval_start) & missing(interval_end)){
      interval_start <- interval_end <- seq(1, length(histogram_data), 1)
    } else if (missing(interval_start)){
      interval_start <- interval_end
    } else if (missing(interval_end)){
      interval_end <- interval_start
    }
  }

  if (!is.double(histogram_data)) {
    histogram_data <- as.double(histogram_data)
  }
  if (!is.integer(interval_start)){
    interval_start <- as.integer(interval_start)
  }
  if (!is.integer(interval_end)){
    interval_end <- as.integer(interval_end)
  }
  if(!is.character(chr)){
    chr <- as.character(chr)
  }

  # Region ID
  if( missing(region_id) & length(histogram_data) > 0 ){
    region_id <- paste0(chr, ":", interval_start[1], "-", interval_end[length(histogram_data)], ":", strand)
  }
  if(!is.character(region_id)){
    region_id <- as.character(region_id)
  }

  # Validate and return object
  return(
    validate_GenomicHistogram(
      new_GenomicHistogram(
        histogram_data = histogram_data,
        interval_start = interval_start,
        interval_end = interval_end,
        chr = chr,
        strand = strand,
        region_id = region_id
      )
    )
  )

}

#' @export
`[.GenomicHistogram` = function(x, i){
  new_GenomicHistogram(
    histogram_data = x$histogram_data[i],
    interval_start = x$interval_start[i],
    interval_end = x$interval_end[i],
    region_id = x$region_id,
    chr = x$chr,
    strand = x$strand)
}


#' @export
reassign_region_id.GenomicHistogram = function(histogram_obj, region_id){

  stopifnot(inherits(histogram_obj, "Histogram"))

  # Creating a region id
  if(missing(region_id) & length(histogram_obj$histogram_data) > 0){
    region_id <- paste0(histogram_obj$chr, ":", histogram_obj$interval_start[1], "-", histogram_obj$interval_end[length(histogram_obj$histogramm_data)], ":", histogram_obj$strand)
  }
  if(!is.character(region_id)){
    region_id <- as.character(region_id)
  }

  histogram_obj$region_id <- region_id

  return(histogram_obj)
}
