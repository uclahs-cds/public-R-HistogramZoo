#' Take a vector of values and get the histogram for integer breaks
#'
#' @param x a vector of observations
#' @param histogram_bin_width Size of histogram bin, assuming a base 1 system
#'
#' @return a Histogram object representing counts of the data
#'
#' @export
#'
#' @examples \dontrun{ observations_to_histogram(sample(1:10, 100, replace = T)) }
observations_to_histogram <- function(x, histogram_bin_width = 1) {

  # Error checking
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(histogram_bin_width))

  # Creating intervals
  a <- floor(min(x))
  b <- ceiling(max(x))

  breaks <- seq(a, b, by = histogram_bin_width)
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
  break_start <- breaks[1:(length(breaks) - 1)]
  break_end <- breaks[2:length(breaks)]

  # Generating a new Histogram
  return(
    new_Histogram(
      histogram_data = as.numeric(hist_data),
      interval_start = break_start,
      interval_end = break_end,
      bin_width = histogram_bin_width,
      region_id = paste0(break_start[1], "-", break_end[length(break_end)])
    )
  )
}
