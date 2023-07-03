#' Return a string representation of an object
#'
#' @param x an object
#'
#' @return a character representation of `x`
dput_str <- function(x) {
  return(
    paste0(utils::capture.output(dput(x)), collapse = " ")
  )
}


#' Correcting for Jaccard/Intersection
#'
#' @param met metric
#' @param value metric fitted value
#'
#' @return if metric is `jaccard` or `intersection`, return 1 - value
correct_fitted_value <- function(met, value){
  if (! is.finite(value)) {
    value <- NA
  } else if (met %in% c("jaccard", "intersection")) {
    return(1 - value)
  } else {
    return(value)
  }
}

#' Convert a vector of points into a data.frame of start/end points representing
#' disjoint intervals
#'
#' @param p integer, a vector of points to be broken up into intervals
#' @param right logical, whether the points should represent interval starts
#' (FALSE) or interval ends (TRUE), default TRUE
#'
#' @return A data.frame with column: start and end representing the start and
#' end points of the intervals
#'
#' @examples \dontrun{
#' index_to_start_end(c(1,5,10))
#' index_to_start_end(c(1,5,10), right = FALSE)
#' }
#'
#' @export
index_to_start_end <- function(p, right = TRUE) {

  # Error checking
  stopifnot(is.numeric(p))
  stopifnot(is.logical(right))

  # Generating segments
  n <- length(p)
  if(n <= 1) {
    stop("Need more than 1 point to compute start/end")
  } else if(n == 2) {
    return_list <- list(
      "start" = p[1],
      "end" = p[2]
    )
  } else if (right) {
    return_list <- list(
      "start" = c(p[1], p[2:(n-1)]+1),
      "end" = p[2:n]
    )
  } else {
    return_list <- list(
      "start" = p[1:(n-1)],
      "end" = c(p[2:(n-1)]-1, p[n])
    )
  }

  return( as.data.frame(return_list) )
}

#' Checks if a numeric vector can be used as an integer vector
#'
#' @param x A numeric vector
#'
#' @return logical, whether `x` is an integer or integer vector
is_equal_integer <- function(x){
  all(x%%1 == 0)
}


#' Finds the bin_width of Histogram or GenomicHistogram bins
#'
#' @param x a Histogram or GenomicHistogram object
#'
#' @return a numeric vector representing bin_widths of each bin
#'
#' @rdname find_bin_width
#' @export
find_bin_width <- function(x){
  UseMethod('find_bin_width')
}

#' @rdname find_bin_width
#' @exportS3Method find_bin_width Histogram
find_bin_width.Histogram <- function(x){
  (x$interval_start + x$interval_end)/2
}

#' @rdname find_bin_width
#' @exportS3Method find_bin_width GenomicHistogram
find_bin_width.GenomicHistogram <- function(x){
  x$consecutive_end - x$consecutive_start + 1
}

#' Finds the midpoint of Histogram or GenomicHistogram bins
#'
#' @param x a Histogram or GenomicHistogram object
#'
#' @return a numeric vector representing midpoints of each bin
#'
#' @rdname find_midpoint
#' @export
find_midpoint <- function(x){
  UseMethod('find_midpoint')
}

#' @rdname find_midpoint
#' @exportS3Method find_midpoint Histogram
find_midpoint.Histogram <- function(x){
  (x$interval_start + x$interval_end)/2
}

#' @rdname find_midpoint
#' @exportS3Method find_midpoint GenomicHistogram
find_midpoint.GenomicHistogram <- function(x){
  (x$consecutive_start + x$consecutive_end)/2
}

#' Loads base config file and overwrites with an optional additional file
#'
#' @param file optional file path to user config script
#'
#' @return list of yaml config
load.config <- function(file = NULL) {
  config.path <- fs::path_package('HistogramZoo', 'config.yaml');
  base.yaml <- yaml::read_yaml(file = config.path);
  if (!is.null(file)) {
    user.yaml <- yaml::read_yaml(file = file);
    base.yaml <- modifyList(base.yaml, user.yaml);
    }

  hz.home <- Sys.getenv('HZ_HOME');
  if (hz.home != '') {
    base.yaml$root.path <- hz.home
    }

  return(base.yaml);
  }

#' Intervals [a, b], [c,d]
#' Need a <= b, c <= d
#' @return TRUE if [a,b] and [c,d] overlap
int_overlap <- function(a,b,c,d) {
  stopifnot(a <= b & c <= d)
  c <= b && a <= d
}
