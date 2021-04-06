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
