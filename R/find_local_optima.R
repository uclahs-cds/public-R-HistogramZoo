

#' find_local_optima returns the local optima of histograms
#'
#' @param x numeric vector representing the density of a histogram
#' @param threshold threshold for local optima, i.e. a point can only be considered a local optima if it differs from its neighbour optima by greater than the permitted threshold, default 0
#' @param flat_endpoints in regions of flat density, whether to return the endpoints or the midpoints
#'
#' @return a list of two vectors
#' \describe{
#'     \item{min_ind}{indices of local minima}
#'     \item{max_ind}{indices of local maxima}
#' }
#'
#' @export
#'
#' @examples \dontrun{
#' x <- c(1,2,3,2,1,2,3,2,1)
#' find_local_optima(x)
#' }
find_local_optima <- function(x, threshold = 0, flat_endpoints = T){
  UseMethod('find_local_optima')
}

#' @exportS3Method find_local_optima Histogram
find_local_optima.Histogram <- function(x, threshold = 0, flat_endpoints = T){
  find_local_optima.numeric(x$histogram_data, threshold = threshold, flat_endpoints = flat_endpoints)
}

#' @exportS3Method find_local_optima numeric
find_local_optima.numeric <- function(x, threshold = 0, flat_endpoints = T){

  # Error checking
  stopifnot(length(x) > 1)
  if(!is_equal_integer(threshold) | !(threshold >= 0) | length(threshold) != 1){
    stop("optima_threshold must be functional as a integer greater than or equal to 0 and of length 1")
  }
  if(!is.logical(flat_endpoints) | length(flat_endpoints) != 1){
    stop("flat_endpoints has to be a logical of length 1")
  }

  # Extracting the changing values in x
  unflat_x <- rle(x)

  # Edge case: If x is a flat block
  if(length(unflat_x$lengths) <= 1) {
    x_min <- c(1, length(x))
    if(!flat_endpoints){
      x_min <- round(mean(x_min))
    }
    return(
      list('min_ind' = x_min, 'max_ind' = NULL)
    )
  }

  # Assuming a stretch of points
  unflat_idx <- cumsum(unflat_x$lengths)
  length_unflat <- length(unflat_idx)
  x_change <- diff(unflat_x$values)

  # Identifying local optima
  x_sign <- sign(x_change)
  x_min <- which(diff(x_sign) == 2)+1
  x_max <- which(diff(x_sign) == -2)+1

  # Filter by threshold
  optima <- sort(c(x_min, x_max))
  length_optima <- length(optima)
  if(length_optima > 1 & threshold > 0){
    optima_vals <- unflat_x$values[optima]
    optima_thresh <- (abs(diff(optima_vals)) > threshold)
    optima_rle <- rle(optima_thresh)
    # if even numbers of consecutive optima points that fail the threshold, deleting them will result in consecutive mins/maxs
    # therefore, change the last consecutive optima pair to TRUE to ensure that consecutive optima points
    # will always be of an odd length
    issue_points <- optima_rle$lengths %% 2 == 0 & !optima_rle$value
    optima_thresh[cumsum(optima_rle$lengths)[issue_points]] <- TRUE
    optima <- optima[c(optima_thresh, TRUE) & c(TRUE, optima_thresh)]
    length_optima <- length(optima)
    x_min <- x_min[x_min %in% optima]
    x_max <- x_max[x_max %in% optima]
  }

  # Filter by (local) threshold
  # x_change_abs <- (abs(x_change) > threshold)
  # x_min <- x_min[x_change_abs[x_min-1] & x_change_abs[x_min]]
  # x_max <- x_max[x_change_abs[x_max-1] & x_change_abs[x_max]]

  # If flat endpoints, take both ends of the flat segments as indices
  if(flat_endpoints){
    x_min_adj <- unique(c(unflat_idx[x_min], unflat_idx[x_min - 1] + 1))
    x_max_adj <- unique(c(unflat_idx[x_max - 1] + 1, unflat_idx[x_max]))
  } else { # Take the midpoint
    x_min_adj <- rowMeans(cbind(unflat_idx[x_min], unflat_idx[x_min - 1] + 1))
    x_min_adj <- round(x_min_adj)
    x_max_adj <- rowMeans(cbind(unflat_idx[x_max], unflat_idx[x_max - 1] + 1))
    x_max_adj <- round(x_max_adj)
  }

  # Start & end points
  start_point <- unique(c(1, unflat_idx[1]))
  end_point <- unique(c(unflat_idx[length_unflat-1]+1, length(x)))
  if(!flat_endpoints) { # Take the midpoint
    start_point <- round(mean(start_point))
    end_point <- round(mean(end_point))
  }

  # Assign start points to min/max
  if ( (length_optima > 0 & optima[1] %in% x_max) |
       (length_optima == 0 & unflat_x$values[2] > unflat_x$values[1]) ){
    x_min_adj <- c(start_point, x_min_adj)
  } else {
    x_max_adj <- c(start_point, x_max_adj)
  }

  # Assign end points to min/max
  if ( (length_optima > 0 & optima[max(length_optima, 1)] %in% x_max) |
       (length_optima == 0 & unflat_x$values[length_unflat] < unflat_x$values[length_unflat - 1]) ){
    x_min_adj <- c(x_min_adj, end_point)
  } else {
    x_max_adj <- c(x_max_adj, end_point)
  }

  # Return results
  return(
    list('min_ind' = x_min_adj, 'max_ind' = x_max_adj)
  )

}
