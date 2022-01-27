#' Return a default value if a variable is null
#'
#' @param a The possibly null value
#' @param b The default value to return if `a` is null
#'
#' @return `a` if `a` is non-null, `b` if `a` is null
#'
#' @examples
#' x <- list(a = 1)
#' y <- x$b %||% 2 # y = 2
#' z <- x$a %||% 10 # z = 1
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}


#' Returns the group matches from a regular expression on a vector
#' @param x the vector we want to match on
#' @param pattern regular expression with groups
str.match <- function(x, pattern) {
  regmatches(x, regexec(pattern, x));
}

#' Return a string representation of an object
dput.str <- function(x) {
  paste0(capture.output(dput(x)), collapse = " ")
}

#' Convert a vector of points into a list of start/end points
#'
#' @param p a vector of indices
#'
#' @return A list with keys: start and end representing the indices
#' @export
#'
#' @examples
#' points.to.start.end(c(1,5,10))
points.to.start.end <- function(p) {
  n = length(p)
  if(n <= 1) {
    stop("Need more than 1 point to compute start/end")
  }
  return.list = list(
    start = p[1:(n - 1)]
  )
  if(n == 2) {
    return.list$end = p[2]
  } else {
    return.list$end = c(p[2:(n - 1)] - 1,  p[length(p)])
  }

  return.list
}
