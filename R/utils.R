#' Returns the group matches from a regular expression on a vector
#' @param x the vector we want to match on
#' @param pattern regular expression with groups
str.match <- function(x, pattern) {
  regmatches(x, regexec(pattern, x));
}

# Return a string representation of an object
dput.str <- function(x) {
  paste0(utils::capture.output(dput(x)), collapse = " ")
}

# Reset rownames
reset.rownames <- function(x) {
  rownames(x) <- NULL
  x
}

#' Convert a vector of points into a list of start/end points
#'
#' @param p a vector of indices
#'
#' @return A data.frame with column: start and end representing the indices
#'
#' @examples \dontrun{
#' index_to_start_end(c(1,5,10))
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
