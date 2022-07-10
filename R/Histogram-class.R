
# constructor
new_Histogram = function(x, interval_start = NULL, interval_end = NULL){
  
  # Checking types
  stopifnot(is.double(x))
  stopifnot(is.integer(interval_start))
  stopifnot(is.integer(interval_end))
  
  # Creating object
  structure(
    x, 
    interval_start = interval_start,
    interval_end = interval_end,
    class = "Histogram"
  )
}

# validator
# TODO: What is the minimum length for a histogram?
# TODO: Do all intervals need to be same length?
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
  
  if(len(values) > 1){
    # 3. Start/End have to be in order
    if(!all(diff(interval_start) > 0 & diff(interval_end) > 0)){
      stop("Intervals must be ordered and nonoverlapping.", call. = FALSE)
    }
    
    # 4. Intervals need to be non-overlapping
    if(any(interval_start[2:(histogram_length)] < interval_end[1:(histogram_length-1)])){
      stop("Intervals must be ordered and nonoverlapping.", call. = FALSE)
    }
  }

  x
}

# helper
Histogram = function(x = double(), interval_start = 1:length(x), interval_end = 1:length(x)){
  
  # Coercing values to the right thing
  x = as.double(x)
  interval_start = as.integer(interval_start)
  interval_end = as.integer(interval_end)
  
  # Validate and return object
  validate_Histogram(new_Histogram(x = x, interval_start = interval_start, interval_end = interval_end))
  
}