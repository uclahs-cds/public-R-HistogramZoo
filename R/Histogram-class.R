
# constructor
new_Histogram = function(x, interval_start = NULL, interval_end = NULL, region_id = NULL, class = character(), ...){

  # Checking types
  stopifnot(is.double(x))
  stopifnot(is.integer(interval_start))
  stopifnot(is.integer(interval_end))
  stopifnot(is.character(region_id))

  # Creating object
  structure(
    x,
    interval_start = interval_start,
    interval_end = interval_end,
    region_id = region_id,
    ...,
    class = c(class, "Histogram")
  )
}

# validator
validate_Histogram = function(x){

  # Attributes
  values = unclass(x)
  histogram_length = length(values)
  interval_start = attr(x, "interval_start")
  interval_end = attr(x, "interval_end")
  region_id = attr(x, "region_id")

  # Validate
  # 0. Everything has to be the same length
  if(length(interval_start) != histogram_length | length(interval_end) != histogram_length) {
    stop("Interval start/end have to be the same length as x", call. = FALSE)
  }

  # 1. Values always have to be greater than 0
  if(!all(!is.na(values) & values >= 0)){
    stop("`x` must be non-missing and nonnegative", call. = FALSE)
  }

  # 2. Start always has to be be less than or equal to end
  if(!all(interval_start <= interval_end)){
    stop("Interval start must be less than interval end for a valid interval", call. = FALSE)
  }

  if(length(values) > 1){
    # 3. Start/End have to be in order
    if(!all(diff(interval_start) > 0 & diff(interval_end) > 0)){
      stop("Intervals must be ordered and nonoverlapping.", call. = FALSE)
    }

    # 4. Intervals need to be non-overlapping
    if(any(interval_start[2:(histogram_length)] <= interval_end[1:(histogram_length-1)])){
      stop("Intervals must be ordered and nonoverlapping.", call. = FALSE)
    }
  }

  # 5. region_id is of length 1
  if(length(region_id) != 1){
    stop("region_id must have length 1.", call. = FALSE)
  }

  x
}

# helper
#' Generates an S3 `Histogram` object
#'
#' @param x vector of counts/density
#' @param interval_start integer vector representing the starts of intervals
#' @param interval_end integer vector representing the ends of intervals
#' @param region_id unique character id to denote region
#'
#' @return A Histogram object
#' @export
#'
#' @examples
#' x = Histogram(x = runif(10), interval_start = 1:10, interval_end = 1:10)
Histogram = function(x = double(), interval_start = integer(), interval_end = integer(), region_id = character()){

  # Coercing values to the right thing
  if(length(x) > 0){
    if( missing(interval_start) & missing(interval_end)){
      interval_start = interval_end = seq(1, length(x), 1)
    } else if (missing(interval_start)){
      interval_start = interval_end
    } else if (missing(interval_end)){
      interval_end = interval_start
    }
    if( missing(region_id) ){
      region_id = paste0(interval_start[1], "-", interval_end[length(x)])
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
  if(!is.character(region_id)){
    region_id = as.character(region_id)
  }

  # Validate and return object
  validate_Histogram(new_Histogram(x = x, interval_start = interval_start, interval_end = interval_end, region_id = region_id))

}

print.Histogram = function(x){

  region_id = attr(x, "region_id")
  interval_start = attr(x, "interval_start")
  interval_end = attr(x, "interval_end")
  models = attr(x, "models")

  # extracting x data
  # TODO: Check with Stefan if there's a more efficient way
  x = unclass(x)
  attributes(x) <- NULL

  # Base case
  cat("Region: ", region_id, "\n")

  # Indices
  if(length(x) > 10){

    # Intervals
    intervals_begin = ifelse(
      interval_start[1:5] == interval_end[1:5],
      interval_start[1:5],
      paste0(interval_start[1:5], "-", interval_end[1:5]))

    intervals_finish = ifelse(
      tail(interval_start, 5) == tail(interval_end, 5),
      tail(interval_start, 5),
      paste0(tail(interval_start, 5), "-", tail(interval_end, 5))
    )
    x_start = as.character(formatC(x[1:5], digits = 2))
    x_end = as.character(formatC(tail(x, 5), digits = 2))

    # Spacing
    spacing = max(nchar(c(intervals_begin, intervals_finish, x_start, x_end)))
    intervals_begin = stringr::str_pad(intervals_begin, width = spacing)
    intervals_finish = stringr::str_pad(intervals_finish, width = spacing)
    x_start = stringr::str_pad(x_start, width = spacing)
    x_end = stringr::str_pad(x_end, width = spacing)

    # Printing
    cat(intervals_begin, "...", intervals_finish, "\n")
    cat(x_start, "...", x_end, "\n")

  } else {

    # Intervals
    intervals = ifelse(interval_start == interval_end, interval_start, paste0(interval_start, "-", interval_end))
    intervals = as.character(intervals)

    # Spacing
    x = as.character(formatC(x, digits = 2))
    spacing = max(nchar(c(intervals, x)))
    intervals = stringr::str_pad(intervals, width = spacing)
    x = stringr::str_pad(x, width = spacing)

    # Printing
    cat(intervals, "\n")
    cat(x, "\n")
  }

  if(!is.null(models)){
    cat("SegmentAndFit: \n")
    cat("Number of Peaks: ", length(models), "\n")
    # TODO: Add stuff about metrics, distributions, parameters, etc
  }

  invisible(x)
}

`[.Histogram` = function(x, i){
  new_Histogram(
    NextMethod(),
    interval_start = attr(x, "interval_start")[i],
    interval_end = attr(x, "interval_end")[i],
    region_id = attr(x, "region_id"))
}

ReassignRegionID = function(x, region_id){
  UseMethod('ReassignRegionID', x)
}

ReassignRegionID.Histogram = function(x, region_id){

  stopifnot(inherits(x, "Histogram"))

  # Creating a region id
  if(missing(region_id)){
    region_id = paste0(attr(x, "interval_start")[1], "-", attr(x, "interval_end")[length(x)])
  }
  if(!is.character(region_id)){
    region_id = as.character(region_id)
  }

  attr(x, "region_id") <- region_id

  return(x)
}
