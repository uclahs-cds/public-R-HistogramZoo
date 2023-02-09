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
  return(
    if (met %in% c("jaccard", "intersection")) (1 - value) else value
  )
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


#' Generate observations approximating a histogram
#'
#' @param x A numeric vector representing the density of a histogram
#' @param rescale_precision numeric, to aid in rescaling fractional densities
#' @return A numeric vector representing approximate observations of that histogram where
#' observations are integer indices of the histogram bin
#' NOTE: This also works on non-integer densities, values are rounded, use with caution
histogram_to_approximate_observations <- function(x, rescale_precision=1){
  return(
    rep(1:length(x), x*rescale_precision)
  )
}
