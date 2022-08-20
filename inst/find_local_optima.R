#' Find changepoints in a vector with uniform stretches of values
#' @param x A numeric vector
find_stepfunction_chgpts <- function(x){
  change_points <- which(diff(x) != 0)
  change_points_plus <- change_points+1
  keep <- (x[change_points] < x[change_points_plus])
  return( c(change_points[keep], change_points_plus[!keep]) )
}

#' Finds the local minima m and maxima M such that
#' m_1 < M_1 < m_2 < M_2 < ... < M_{K - 1} < m_{k}
#'
#' @param x a numeric vector, representing the density of a histogram
#' @param threshold numeric, minimum distance between local optima
#'
#' @export
find_local_optima <- function(x, threshold = 0) {

  # Error checking
  x <- as.numeric(x)
  stopifnot(length(x) > 1)

  # Get the first non-equal index
  n <- length(x)
  min_ind <- NULL
  max_ind <- NULL
  #init = x[1]
  # If x = c(1,1,1,1,2,3) then j = 4, the first index before the run ends

  # Trim the left and right side
  rle_diff <- rle(diff(x))
  # Return no min/max if all equal
  if(length(rle_diff$lengths) <= 1) {
    return(
      list(
        'min_ind' = 1,
        'max_ind' = NULL
        )
      )
  }

  start_index <- 1
  end_index <- n
  if(rle_diff$values[1] == 0) {
    start_index <- rle_diff$lengths[1] + 1
  }
  if(length(rle_diff$lengths) > 2 && utils::tail(rle_diff$values, n = 1) == 0) {
    # Remove the last equal values
    end_index <- n - utils::tail(rle_diff$lengths, n = 1)
  }
  x_trim <- x[start_index:end_index]
  n_trim <- length(x_trim)
  # Keep track if we last appended a min/max
  min_appended <- NULL
  # Check the first point for min/max
  if(x_trim[1] < x_trim[2]) {
    min_ind <- 1
    min_appended <- TRUE
  } else if (x_trim[1] > x_trim[2]) {
    max_ind <- 1
    min_appended <- FALSE
  }

  # Do the middle segment
  if(n_trim > 3) {
    for(i in seq(2, n_trim - 1)) {
      # min_appended ensures that we alternate minima and maxima
      left_diff <- x_trim[i] - x_trim[i - 1]
      right_diff <- x_trim[i] - x_trim[i + 1]

      if (!min_appended &&
          ((left_diff < -threshold && right_diff <= 0) ||
           (left_diff <= 0 && right_diff < -threshold))) {
        min_ind <- c(min_ind, i)
        min_appended <- TRUE
      } else if(min_appended &&
                ((left_diff > threshold && right_diff >= 0) ||
                (left_diff >= 0 && right_diff > threshold))) {
        max_ind <- c(max_ind, i)
        min_appended <- FALSE
      }
    }
  }

  if(x_trim[n_trim - 1] > x_trim[n_trim]) {
    min_ind <- c(min_ind, n_trim)
    min_appended <- TRUE
  } else if (x_trim[n_trim - 1] < x_trim[n_trim]) {
    max_ind <- c(max_ind, n_trim)
    min_appended <- FALSE
  }

  return(
    list(
      'min_ind' = min_ind + (start_index - 1),
      'max_ind' = max_ind + (start_index - 1)
    )
  )

}
