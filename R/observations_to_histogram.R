#' Take a vector of values and get the histogram for integer breaks
#'
#' @param x A vector of observations
#' @param histogram_bin_width Size of histogram bin, assuming a base 1 system
#'
#' @return A Histogram object representing counts of the data
#'
#' @export
#'
#' @examples \dontrun{ observations_to_histogram(sample(1:10, 100, replace = T)) }
observations_to_histogram = function(x, histogram_bin_width = 1) {

  # Error checking
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(histogram_bin_width))

  # Creating intervals
  a <- floor(min(x))
  b <- ceiling(max(x))
  breaks <- seq(a - 1, b, by = histogram_bin_width)
  breaks <- unique(c(breaks, b))

  # Histogram data
  hist_data <- table(
    cut(x,
        breaks = breaks,
        right = TRUE,
        include.lowest = TRUE
        )
  )

  # Generating endpoints
  breaks <- c(a, breaks[breaks != (a-1)])
  endpoints <- index_to_start_end(breaks, right = TRUE)

  # Generating a new Histogram
  return(
    new_Histogram(
      histogram_data = as.numeric(hist_data),
      interval_start = as.integer(endpoints$start),
      interval_end = as.integer(endpoints$end),
      region_id = paste0(endpoints$start[1], "-", endpoints$end[length(hist_data)])
    )
  )
}
