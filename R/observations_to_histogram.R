#' Take a vector of values and get the histogram for integer breaks
#'
#' @param x A vector of observations
#' @param histogram_bin_width Maximum size of histogram bin, assuming a base 1 system
#' @param right Logical, indicating if histogram bins should be closed on the right. If FALSE, bins will be closed on the left.
#' @param include_endpoint Whether the boundary value (highest or lowest) should be counted in the bins.
#'
#' @return A Histogram object representing counts of the data
#'
#' @export
#' 
#' @example \dontrun{
#' observations_to_histogram(sample(1:10, 100, replace = T))
#' }
observations_to_histogram = function(x, histogram_bin_width = 1, right = TRUE, include_endpoint = TRUE) {
  
  # Error checking
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(histogram_bin_width))
  stopifnot(is.logical(right))
  stopifnot(is.logical(include_endpoint))
  
  # Creating intervals
  a = floor(min(x))
  b = ceiling(max(x))
  histogram_range = b-a+1
  length_out = ceiling(histogram_range/histogram_bin_width)+1
  breaks = seq(a, b, length.out = length_out)
  breaks = round(breaks)
  
  # Histogram data
  hist_data = table(
    cut(x, 
        breaks = breaks, 
        right = right, 
        include.lowest = include_endpoint)
  )
  
  # Generating endpoints
  if(right){
    endpoints = index_to_start_end(breaks, right = TRUE)
  } else{
    endpoints = index_to_start_end(breaks, right = FALSE)
  }
  
  # Generating a new Histogram
  return(
    new_Histogram(
      histogram_data = hist_data, 
      interval_start = endpoints$start, 
      interval_end = endpoints$end, 
      region_id = NULL
    )
  )
}