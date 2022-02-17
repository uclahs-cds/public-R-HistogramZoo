#' Returns the group matches from a regular expression on a vector
#' @param x the vector we want to match on
#' @param pattern regular expression with groups
str.match <- function(x, pattern) {
  regmatches(x, regexec(pattern, x));
}

# Return a string representation of an object
dput.str <- function(x) {
  paste0(capture.output(dput(x)), collapse = " ")
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
#' index.to.start.end(c(1,5,10))
#' }
#'
#' @export
index.to.start.end <- function(p) {
  n = length(p)
  if(n <= 1) {
    stop("Need more than 1 point to compute start/end")
  }
  return.list = list(
    start = p[1:(n - 1)]
  )
  if(n == 2) {
    return.list = list(start = p[1])
  } else {
    return.list = list(start = c(p[1],  p[2:(length(p)-1)]+1))
  }
  return.list$end = c(p[2:length(p)])

  as.data.frame(return.list)
}
