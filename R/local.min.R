# From internal function ggpmisc:::find_peaks

#' Find local maxima or global maximum (peaks)
#'
#' This method finds peaks (local maxima) in a spectrum, using a user selectable
#' span and size threshold relative to the tallest peak (global maximum). This
#' is a wrapper built on top of function peaks from package 'splus2R'.
#'
#' @param x numeric array
#' @param ignore_threshold numeric value between 0.0 and 1.0 indicating the size
#'   threshold below which peaks will be ignored.
#' @param span a peak is defined as an element in a sequence which is greater
#'   than all other elements within a window of width span centered at that
#'   element. The default value is 3, meaning that a peak is bigger than both of
#'   its neighbors. Default: 3.
#' @param strict logical flag: if TRUE, an element must be strictly greater than
#'   all other values in its window to be considered a peak. Default: TRUE.
#' @param na.rm logical indicating whether \code{NA} values should be stripped
#'   before searching for peaks.
#'
#' @return A vector of logical values. Values that are TRUE
#'   correspond to local peaks in the data.
#'
#' @details This function is a wrapper built on function
#'   \code{\link[splus2R]{peaks}} from \pkg{splus2R} and handles non-finite
#'   (including NA) values differently than \code{peaks}, instead of giving an
#'   error they are replaced with the smallest finite value in \code{x}.
#'
#' @note The default for parameter \code{strict} is \code{TRUE} in functions
#'   \code{peaks()} and \code{find_peaks()}, while the default is \code{FALSE}
#'   in \code{stat_peaks()} and in \code{stat_valleys()}.
#'
#' @seealso \code{\link[splus2R]{peaks}}
#'
#' @family peaks and valleys functions
#'
#' @keywords internal
#'
find_peaks <-
  function(x,
           ignore_threshold = 0,
           span = 3,
           strict = TRUE,
           na.rm = FALSE) {
    # find peaks
    if(is.null(span)) {
      pks <- x == max(x, na.rm = na.rm)
      if (strict && sum(pks) != 1L) {
        pks <- logical(length(x)) # all FALSE
      }
    } else {
      pks <- splus2R::peaks(x = x, span = span, strict = strict)
    }
    # apply threshold to found peaks
    if (abs(ignore_threshold) < 1e-5) {
      pks
    } else {
      range_x <- range(x, na.rm = na.rm, finite = TRUE)
      min_x <- range_x[1]
      max_x <- range_x[2]
      x <- ifelse(!is.finite(x), min_x, x)
      # this can cater for the case when max_x < 0, as with logs
      delta <- max_x - min_x
      top_flag <- ignore_threshold > 0.0
      scaled_threshold <- delta * abs(ignore_threshold)
      if (top_flag) {
        ifelse(x - min_x > scaled_threshold, pks , FALSE)
      } else {
        ifelse(max_x - x > scaled_threshold, pks , FALSE)
      }
    }
  }
