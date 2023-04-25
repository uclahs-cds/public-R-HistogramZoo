#' Find changepoints in a vector with uniform stretches of values
#' @param x A numeric vector or Histogram object
#'
#' @return A numeric vector representing change point indices
#'
#' @rdname find_change_points
#' @export
find_change_points <- function(x){
  UseMethod("find_change_points")
}

#' @rdname find_change_points
#' @exportS3Method find_change_points Histogram
find_change_points.Histogram <- function(x){
  find_change_points.numeric(x$histogram_data)
}

#' @rdname find_change_points
#' @exportS3Method find_change_points numeric
find_change_points.numeric <- function(x){
  change_points <- which(diff(x) != 0)
  change_points_plus <- change_points+1
  keep <- (x[change_points] < x[change_points_plus])
  return( c(change_points[keep], change_points_plus[!keep]) )
}
