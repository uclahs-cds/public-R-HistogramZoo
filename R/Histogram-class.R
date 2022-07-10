
# constructor
new_Histogram = function(x, interval_start = NULL, interval_end = NULL, class = character(), ...){

  # Checking types
  stopifnot(is.double(x))
  stopifnot(is.integer(interval_start))
  stopifnot(is.integer(interval_end))

  # Creating object
  structure(
    x,
    interval_start = interval_start,
    interval_end = interval_end,
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

  # Validate
  # 0. Everything has to be the same length
  if(length(interval_start) != histogram_length | length(interval_end) != histogram_length) {
    stop("Interval start/end have to be the same length as x", call. = FALSE)
  }

  # 1. Values always have to be greater than 0
  if(!all(values > 0 & values > 0)){
    stop("`x` must be non-missing and greater than 0", call. = FALSE)
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

  x
}

# helper
#' Generates an S3 `Histogram` object
#'
#' @param x vector of counts/density
#' @param interval_start integer vector representing the starts of intervals
#' @param interval_end integer vector representing the ends of intervals
#'
#' @return A Histogram object
#' @export
#'
#' @examples
#' x = Histogram(x = runif(10), interval_start = 1:10, interval_end = 1:10)
Histogram = function(x = double(), interval_start = integer(), interval_end = integer()){

  # Coercing values to the right thing
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

  # Validate and return object
  validate_Histogram(new_Histogram(x = x, interval_start = interval_start, interval_end = interval_end))

}
