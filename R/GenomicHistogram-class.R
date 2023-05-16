
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
#' @param intron_start integer vector representing the starts of introns
#' @param intron_end integer vector representing the ends of introns
#' @param consecutive_start start of intervals, integer vector representing an intronless set of bins - defaults to starting at 1
#' @param consecutive_end end of intervals, integer vector representing an intronless set of bins - defaults to starting at `bin_width`
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
    intron_start = NULL,
    intron_end = NULL,
    consecutive_start = NULL,
    consecutive_end = NULL
  ){

  # Checking types
  stopifnot(is.character(chr))
  stopifnot(is.character(strand))
  stopifnot(is.numeric(intron_start))
  stopifnot(is.numeric(intron_end))
  stopifnot(is.numeric(consecutive_start))
  stopifnot(is.numeric(consecutive_end))

  # Creating object
  return(
    new_Histogram(
      histogram_data = histogram_data,
      interval_start = interval_start,
      interval_end = interval_end,
      bin_width = bin_width,
      region_id = region_id,
      chr = chr,
      strand = strand,
      intron_start = intron_start,
      intron_end = intron_end,
      consecutive_start = consecutive_start,
      consecutive_end = consecutive_end,
      class = "GenomicHistogram"
    )
  )
}

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
  bin_gr <- IRanges::IRanges(start = interval_start, end = interval_end)
  range_gr <- base::range(bin_gr)

  # Intron attributes
  intron_start <- x$intron_start
  intron_end <- x$intron_end
  intron_length <- length(intron_start)
  intron_gr <- IRanges::IRanges(start = intron_start, end = intron_end)

  # Consecutive start and end
  consecutive_start <- x$consecutive_start
  consecutive_end <- x$consecutive_end

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

  # 5. region_id is of length 1
  if(length(region_id) != 1){
    stop("region_id must have length 1.", call. = FALSE)
  }

  # 6. chr is of length 1
  if( length(chr) > 1 ) {
    stop("chr must have length 1", call. = FALSE)
  }

  # 7. Introns
  if(length(intron_start) != length(intron_end)){
    stop("intron start/end have to be the same length", call. = FALSE)
  }

  if(!all(intron_start <= intron_end)){
    stop("intron_start must be less than or equal to intron_end for a valid intron", call. = FALSE)
  }

  if(intron_length > 1){
    # 8. Start/End have to be in order
    if(!all(diff(intron_start) > 0 & diff(intron_end) > 0)){
      stop("introns must be ordered and nonoverlapping.", call. = FALSE)
    }

    # 9. Introns need to be non-overlapping
    if(any(intron_start[2:(intron_length)] <= intron_end[1:(intron_length-1)])){
      stop("introns must be ordered and nonoverlapping.", call. = FALSE)
    }
  }

  if(intron_length > 0){
    # 10. All estimated introns are actually in the set of introns
    intron_est <- IRanges::setdiff(range_gr, bin_gr)
    if(length(intron_est) > 0){
      ovl <- IRanges::findOverlaps(intron_est, intron_gr, type = "equal")
      if(!all(seq(1, length(intron_est)) %in% S4Vectors::queryHits(ovl))){
        stop("estimated introns not found in intron set", call. = FALSE)
      }
    }

    # 11. Check that every intron is in the range
    ovl <- IRanges::findOverlaps( intron_gr, range_gr, type = "within")
    if(any(!seq(1, intron_length) %in% S4Vectors::queryHits(ovl))){
      stop("introns must overlap histogram range", call. = FALSE)
    }
  }

  # 12. bin_width is positive integer
  if(!(is.integer(bin_width) && bin_width > 0)){
    stop("bin_width must be a positive integer", call. = FALSE)
  }

  # 13. All histogram bins are length bin_width
  bin_split <- split(bin_gr, seq(1, histogram_length))
  bin_width_vec <- sapply(bin_split, function(gr) sum(IRanges::width(IRanges::setdiff(gr, intron_gr))))

  if(!all(bin_width_vec[1:(histogram_length-1)] == bin_width)){
    stop(
      'Incorrect bin_width after intron correction:
       If GenomicHistogram is of greater than length 1, all bins subtracting introns
       with the exception of the last bin must be equal in length to bin_width.
       For a GenomicHistogram of length 1, the bin length of the only bin must be equal to bin_width',
      call. = FALSE
    )
  }

  # 14. Checking basic properties of consecutive bin width
  if(!all(consecutive_start <= consecutive_end)){
    stop("consecutive_start must be less than or equal to consecutive_end for a valid interval", call. = FALSE)
  }

  if(histogram_length > 1){
    if(!all(diff(consecutive_start) > 0 & diff(consecutive_end) > 0)){
      stop("consecutive intervals must be ordered and nonoverlapping.", call. = FALSE)
    }

    if(!all(consecutive_start[2:(histogram_length)] == consecutive_end[1:(histogram_length-1)] + 1)){
      stop("consecutive intervals must be consecutive", call. = FALSE)
    }
  }

  # 15. Consecutive bins follow estimated bin width
  consecutive_bin_vec <- consecutive_end - consecutive_start + 1
  if(!all(consecutive_bin_vec == bin_width_vec)){
    stop('consecutive bin width not equal to estimated bin width, missing introns or incorrect consecutive bin specification')
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
#' @param intron_start integer vector representing the starts of introns (optional: to represent intron-spanning histograms on transcripts)
#' @param intron_end integer vector representing the ends of introns (optional: to represent intron-spanning histograms on transcripts)
#' @param consecutive_start start of intervals, integer vector representing an intronless set of bins - defaults to starting at 1
#' @param consecutive_end end of intervals, integer vector representing an intronless set of bins - defaults to starting at `bin_width`
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
    strand = c("*", "+", "-"),
    intron_start = integer(),
    intron_end = integer(),
    consecutive_start = integer(),
    consecutive_end = integer()
){

  # Coercing values to the right thing
  strand <- match.arg(strand)

  # Coercing to the right data type
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
  if(!is.integer(intron_start)){
    intron_start <- as.integer(intron_start)
  }
  if(!is.integer(intron_end)){
    intron_end <- as.integer(intron_end)
  }
  if(!is.character(region_id)){
    region_id <- as.character(region_id)
  }
  if(!is.integer(consecutive_start)){
    consecutive_start <- as.integer(consecutive_start)
  }
  if(!is.integer(consecutive_end)){
    consecutive_end <- as.integer(consecutive_end)
  }

  # Estimating missing parameters
  if(length(histogram_data) > 0){
    # Assigning interval start and end if missing
    if( missing(interval_start) & missing(interval_end)){
      interval_start <- interval_end <- seq(1, length(histogram_data), 1)
    } else if (missing(interval_start)){
      interval_start <- interval_end
    } else if (missing(interval_end)){
      interval_end <- interval_start
    }

    bins <- IRanges::IRanges(start = interval_start, end = interval_end)

    # Adding missing introns
    range_gr <- base::range(bins)
    introns <- IRanges::setdiff(range_gr, bins)
    intron_start <- c( BiocGenerics::start(introns), intron_start)
    intron_end <- c( BiocGenerics::end(introns), intron_end)
    intron_dup <- !duplicated(cbind(intron_start, intron_end))
    intron_start <- intron_start[intron_dup]
    intron_end <- intron_end[intron_dup]

    # Reordering
    intron_start_order <- order(intron_start)
    intron_end_order <- order(intron_end)
    if(!all(intron_start_order == intron_end_order)){
      stop("Something wrong with introns - likely overlap - revise GenomicHistogram")
    }
    intron_start <- intron_start[intron_start_order]
    intron_end <- intron_end[intron_end_order]

    # Estimating bin width if missing
    if(missing(bin_width)){
      intron_gr <- IRanges::IRanges(start = intron_start, end = intron_end)
      bin_width <- as.integer(
        sum(IRanges::width(IRanges::setdiff(bins[1], intron_gr)))
      )
    }

    # Estimating consecutive intervals
    if(missing(consecutive_start) | missing(consecutive_end)){

      consecutive_start <- seq(0, length(histogram_data)-1)*bin_width + 1
      consecutive_end <- seq(1, length(histogram_data))*bin_width

      # Estimating last bin width
      last_bin_intron <- (intron_start >= tail(interval_start, 1) & intron_end <= tail(interval_end, 1))
      last_bin_width <- (tail(interval_end, 1) - tail(interval_start, 1) + 1) - sum(intron_end[last_bin_intron] - intron_start[last_bin_intron] + 1)
      consecutive_end[length(histogram_data)] <- consecutive_start[length(histogram_data)] + last_bin_width - 1
    }

  }

  # Region ID
  if( missing(region_id) & length(histogram_data) > 0 ){
    region_id <- paste0(chr, ":", interval_start[1], "-", interval_end[length(histogram_data)], ":", strand)
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
        strand = strand,
        intron_start = intron_start,
        intron_end = intron_end,
        consecutive_start = consecutive_start,
        consecutive_end = consecutive_end
      )
    )
  )

}

#' Formats print intervals with introns and one-bp bins
format_print_intervals <- function(
  interval_start,
  interval_end,
  bin_size_one,
  intron_bins
){
  x <- vector(mode="character", length=length(interval_start))
  x[bin_size_one] <- interval_start[bin_size_one]
  x[!bin_size_one & intron_bins] <- paste0(
    interval_start[!bin_size_one & intron_bins],
    "><",
    interval_end[!bin_size_one & intron_bins]
  )
  x[!bin_size_one & !intron_bins] <- paste0(
    interval_start[!bin_size_one & !intron_bins],
    "-",
    interval_end[!bin_size_one & !intron_bins]
  )
  return(x)
}

#' @export
print.GenomicHistogram <- function(x, ...){

  histogram_data <- x$histogram_data
  region_id <- x$region_id
  interval_start <- x$interval_start
  interval_end <- x$interval_end
  intron_start <- x$intron_start
  intron_end <- x$intron_end

  # Base case
  cat("Region: ", region_id, "\n")

  # Identifying introns within bins
  intron_bins <- sapply(
    seq_along(interval_start),
    function(i) any(intron_start >= interval_start[i] & intron_end <= interval_end[i])
  )

  # Representation bin-size 1 - this doesn't mean functionally bin-size 1
  bin_size_one <- interval_start == interval_end

  # Indices
  if(length(histogram_data) > 10){

    # Intervals
    intervals_begin <- format_print_intervals(
      interval_start[1:5],
      interval_end[1:5],
      bin_size_one[1:5],
      intron_bins[1:5]
    )

    intervals_finish <- format_print_intervals(
      utils::tail(interval_start, 5),
      utils::tail(interval_end, 5),
      utils::tail(bin_size_one, 5),
      utils::tail(intron_bins, 5)
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
    intervals <- format_print_intervals(
      interval_start,
      interval_end,
      bin_size_one,
      intron_bins
    )

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
`[.GenomicHistogram` <- function(x, i){

  # Continuity enforcement
  if(length(i) > 1 && !all(diff(i)== 1)){
    stop("Subsetting of Histograms must include a valid continuous set of indices")
  }

  # Keep only introns in the subset
  keep_introns <- x$intron_start > x$interval_start[i][1] & x$intron_end < tail(x$interval_end[i], n=1)

  # Change bin width if necessary
  if(length(i) == 1 && i == length(x)){
    x$bin_width <- x$interval_end[i] - x$interval_start[i] + 1
    if(any(keep_introns)){ # If introns need to be kept
      x$bin_width <- x$bin_width - sum(x$intron_end[keep_introns] - x$intron_start[keep_introns] + 1)
    }
  }

  # Generate new histogram
  new_GenomicHistogram(
    histogram_data = x$histogram_data[i],
    interval_start = x$interval_start[i],
    interval_end = x$interval_end[i],
    region_id = x$region_id,
    bin_width = x$bin_width,
    chr = x$chr,
    strand = x$strand,
    intron_start = x$intron_start[keep_introns],
    intron_end = x$intron_end[keep_introns],
    consecutive_start = x$consecutive_start[i],
    consecutive_end = x$consecutive_end[i]
  )
}


#' @export
reassign_region_id.GenomicHistogram <- function(histogram_obj, region_id){

  stopifnot(inherits(histogram_obj, "Histogram"))

  # Creating a region id
  if(missing(region_id) & length(histogram_obj$histogram_data) > 0){
    region_id <- paste0(histogram_obj$chr, ":", histogram_obj$interval_start[1], "-", histogram_obj$interval_end[length(histogram_obj$histogram_data)], ":", histogram_obj$strand)
  }
  if(!is.character(region_id)){
    region_id <- as.character(region_id)
  }

  histogram_obj$region_id <- region_id

  return(histogram_obj)
}


#' Reset consecutive intervals to starting at 1, bin_width
#'
#' @param x GenomicHistogram object
#'
#' @return a GenomicHistogram object where `consecutive_start` is reset to 1 and `consecutive_end` is reset to bin_width
#'
#' @export
#'
#' @examples \dontrun{
#' x <- GenomicHistogram(c(1,2,3,4,5))
#' x <- x[3:5]
#' reset_consecutive_intervals(x)
#' }
reset_consecutive_intervals <- function(x){

  stopifnot(inherits(x, "GenomicHistogram"))

  new_GenomicHistogram(
    histogram_data = x$histogram_data,
    interval_start = x$interval_start,
    interval_end = x$interval_end,
    region_id = x$region_id,
    bin_width = x$bin_width,
    chr = x$chr,
    strand = x$strand,
    intron_start = x$intron_start,
    intron_end = x$intron_end,
    consecutive_start = x$consecutive_start - head(x$consecutive_start, 1) + 1,
    consecutive_end = x$consecutive_end - head(x$consecutive_start, 1) + 1
  )
}
