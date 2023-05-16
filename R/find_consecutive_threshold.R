#' Returns the indices for consecutive elements of a vector that are greater than a specified threshold
#'
#' @param x A numeric vector or Histogram object
#' @param threshold numeric threshold
#'
#' @return list of coordinates with `start` and `end` coordinates
#'
#' @rdname find_consecutive_threshold
#' @export
#'
#' @examples
#' find_consecutive_threshold(c(0,0,0,1,1,1,0,0,0,1,1,1,0,0))
#' find_consecutive_threshold(c(0,0,1,2,2,0,1,1,1,0,0), threshold = 1)
find_consecutive_threshold <- function(x, threshold = 0){
  UseMethod('find_consecutive_threshold')
}

#' @rdname find_consecutive_threshold
#' @exportS3Method find_consecutive_threshold Histogram
find_consecutive_threshold.Histogram <- function(x, threshold = 0){
  find_consecutive_threshold.numeric(x = x$histogram_data, threshold = threshold)
}

#' @rdname find_consecutive_threshold
#' @exportS3Method find_consecutive_threshold numeric
find_consecutive_threshold.numeric <- function(x, threshold = 0){
  x_thresholded <- rle(x > threshold)
  end_coords <- cumsum(x_thresholded$lengths)
  start_coords <- end_coords - x_thresholded$lengths + 1
  start_coords_thresholded <- start_coords[x_thresholded$values]
  end_coords_thresholded <- end_coords[x_thresholded$values]

  return(list(start = start_coords_thresholded, end = end_coords_thresholded))
}
