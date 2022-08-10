#' Return a string representation of an object
#' @param x an object
dput_str <- function(x) {
  return(
    paste0(utils::capture.output(dput(x)), collapse = " ")
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
#' @param x A numeric vector
is_equal_integer <- function(x){
  all(x%%1 == 0)
}
