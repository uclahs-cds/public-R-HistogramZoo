
# constructor
#' GenomicHistogram constructor
#'
#' @param histogram_data vector of counts/density
#' @param interval_start integer vector representing the starts of intervals
#' @param interval_end integer vector representing the ends of intervals
#' @param region_id character identifier for the region of interest
#' @param bin_width integer width of histogram bins, if missing, estimated from `interval_start` and `interval_end`
#' @param chr chromosome name
#' @param strand strand
#'
#' @return A GenomicHistogram object
new_GenomicHistogram <- function(
    histogram_data = NULL,
    interval_start = NULL,
    interval_end = NULL,
    region_id = NULL,
    bin_width = NULL,
    chr = NULL,
    strand = NULL,
    ...
  ){

  # Checking types
  stopifnot(is.double(histogram_data))
  stopifnot(is.integer(interval_start))
  stopifnot(is.integer(interval_end))
  stopifnot(is.character(region_id))
  stopifnot(is.integer(bin_width))
  stopifnot(is.character(chr))
  stopifnot(is.character(strand))

  # Creating object
  x <- list(
    histogram_data = histogram_data,
    interval_start = interval_start,
    interval_end = interval_end,
    bin_width = bin_width,
    region_id = region_id,
    chr = chr,
    strand = strand,
    ...
  )
  class(x) <- c(class, "GenomicHistogram")

  return(x)
}

# validator
#' Validates GenomicHistogram objects
#'
#' @param x a GenomicHistogram object candidate
#' @return a validated GenomicHistogram object
validate_GenomicHistogram <- function(x){

  # Attributes
  histogram_data <- x$histogram_data
  histogram_length <- length(histogram_data)
  interval_start <- x$interval_start
  interval_end <- x$interval_end
  bin_width <- x$bin_width
  region_id <- x$region_id
  chr <- x$chr

  # Validate
  # 0. Everything has to be the same length, Histogram has to be at least length 1
  if(histogram_length == 0){
    stop("Histogram must be at least length 1", call. = FALSE)
  }

  if(length(interval_start) != histogram_length | length(interval_end) != histogram_length) {
    stop("interval_start/end have to be the same length as histogram_data", call. = FALSE)
  }

  # 1. histogram_data always have to be greater than 0
  if(!all(!is.na(histogram_data) & histogram_data >= 0)){
    stop("histogram_data must be non-missing and nonnegative", call. = FALSE)
  }

  # 2. Start always has to be be less than or equal to end
  if(!all(interval_start <= interval_end)){
    stop("interval_start must be less than or equal to interval_end for a valid interval", call. = FALSE)
  }

  if(histogram_length > 1){
    # 3. Start/End have to be in order
    if(!all(diff(interval_start) > 0 & diff(interval_end) > 0)){
      stop("intervals must be ordered and nonoverlapping.", call. = FALSE)
    }

    # 4. Intervals need to be non-overlapping
    if(any(interval_start[2:(histogram_length)] <= interval_end[1:(histogram_length-1)])){
      stop("intervals must be ordered and nonoverlapping.", call. = FALSE)
    }
  }

  # 5. bin_width is positive integer
  if(!(is.integer(bin_width) && bin_width > 0)){
    stop("bin_width must be a positive integer", call. = FALSE)
  }

  # 6. All histogram bins are length bin_width
  bin_width_vec <- interval_end - interval_start + 1
  if(!all(bin_width_vec[1:(histogram_length-1)] == bin_width)){
    stop("Incorrect bin_width: If Histogram is of greater than length 1, all intervals with the exception of the last interval must be equal in length to bin_width. For a Histogram of length 1, the interval length must be equal to bin_width", call. = FALSE)
  }

  # 7. region_id is of length 1
  if(length(region_id) != 1){
    stop("region_id must have length 1.", call. = FALSE)
  }

  # 8. chr is of length 1
  if( length(chr) > 1 ) {
    stop("chr must have length 1", call. = FALSE)
  }

  return(x)
}

# helper
#' Generates an S3 `GenomicHistogram` object
#'
#' @param histogram_data vector of counts/density
#' @param interval_start integer vector representing the starts of intervals
#' @param interval_end integer vector representing the ends of intervals
#' @param region_id character identifier for the region of interest
#' @param bin_width integer width of histogram bins, if missing, estimated from `interval_start` and `interval_end`
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
    bin_width = integer(),
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
    if(missing(bin_width)){
      bin_width <- as.integer(
        ceiling(max(interval_end - interval_start + 1))
      )
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
  if (!is.integer(bin_width)){
    bin_width <- as.integer(bin_width)
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
        bin_width = bin_width,
        region_id = region_id,
        chr = chr,
        strand = strand
      )
    )
  )

}

#' @export
print.GenomicHistogram = function(x, ...){

  histogram_data <- x$histogram_data
  region_id <- x$region_id
  interval_start <- x$interval_start
  interval_end <- x$interval_end

  # Base case
  cat("Region: ", region_id, "\n")

  # Indices
  if(length(histogram_data) > 10){

    # Intervals
    intervals_begin <- ifelse(
      interval_start[1:5] == interval_end[1:5],
      interval_start[1:5],
      paste0(interval_start[1:5], "-", interval_end[1:5]))

    intervals_finish <- ifelse(
      utils::tail(interval_start, 5) == utils::tail(interval_end, 5),
      utils::tail(interval_start, 5),
      paste0(utils::tail(interval_start, 5), "-", utils::tail(interval_end, 5))
    )
    x_start <- as.character(formatC(histogram_data[1:5], digits = 2))
    x_end <- as.character(formatC(utils::tail(histogram_data, 5), digits = 2))

    # Spacing
    spacing <- max(nchar(c(intervals_begin, intervals_finish, x_start, x_end)))
    intervals_begin <- stringr::str_pad(intervals_begin, width = spacing)
    intervals_finish <- stringr::str_pad(intervals_finish, width = spacing)
    x_start <- stringr::str_pad(x_start, width = spacing)
    x_end <- stringr::str_pad(x_end, width = spacing)

    # Printing
    cat(intervals_begin, "...", intervals_finish, "\n")
    cat(x_start, "...", x_end, "\n")

  } else {
    # Intervals
    intervals <- ifelse(
      interval_start == interval_end,
      interval_start,
      paste0(interval_start, "-", interval_end))

    intervals <- as.character(intervals)

    # Spacing
    x_full <- as.character(formatC(histogram_data, digits = 2))
    spacing <- max(nchar(c(intervals, x_full)))
    intervals <- stringr::str_pad(intervals, width = spacing)
    x_full <- stringr::str_pad(x_full, width = spacing)

    # Printing
    cat(intervals, "\n")
    cat(x_full, "\n")
  }

  return(invisible(x))
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
