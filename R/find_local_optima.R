
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
  min.ind <- NULL
  max.ind <- NULL
  #init = x[1]
  # If x = c(1,1,1,1,2,3) then j = 4, the first index before the run ends

  # Trim the left and right side
  rle_diff <- rle(diff(x))
  # Return no min/max if all equal
  if(length(rle_diff$lengths) <= 1) {
    return(
      list(
        'min.ind' = 1,
        'max.ind' = NULL
        )
      )
  }

  startIndex <- 1
  endIndex <- n
  if(rle_diff$values[1] == 0) {
    startIndex <- rle_diff$lengths[1] + 1
  }
  if(length(rle_diff$lengths) > 2 && utils::tail(rle_diff$values, n = 1) == 0) {
    # Remove the last equal values
    endIndex <- n - utils::tail(rle_diff$lengths, n = 1)
  }
  x.trim <- x[startIndex:endIndex]
  n.trim <- length(x.trim)
  # Keep track if we last appended a min/max
  min.appended <- NULL
  # Check the first point for min/max
  if(x.trim[1] < x.trim[2]) {
    min.ind <- 1
    min.appended <- TRUE
  } else if (x.trim[1] > x.trim[2]) {
    max.ind <- 1
    min.appended <- FALSE
  }

  # Do the middle segment
  if(n.trim > 3) {
    for(i in seq(2, n.trim - 1)) {
      # min.appended ensures that we alternate minima and maxima
      left.diff <- x.trim[i] - x.trim[i - 1]
      right.diff <- x.trim[i] - x.trim[i + 1]

      if (!min.appended &&
          ((left.diff < -threshold && right.diff <= 0) ||
           (left.diff <= 0 && right.diff < -threshold))) {
        min.ind <- c(min.ind, i)
        min.appended <- TRUE
      } else if(min.appended &&
                ((left.diff > threshold && right.diff >= 0) ||
                (left.diff >= 0 && right.diff > threshold))) {
        max.ind <- c(max.ind, i)
        min.appended <- FALSE
      }
    }
  }

  if(x.trim[n.trim - 1] > x.trim[n.trim]) {
    min.ind <- c(min.ind, n.trim)
    min.appended <- TRUE
  } else if (x.trim[n.trim - 1] < x.trim[n.trim]) {
    max.ind <- c(max.ind, n.trim)
    min.appended <- FALSE
  }

  return(
    list(
      'min.ind' = min.ind + (startIndex - 1),
      'max.ind' = max.ind + (startIndex - 1)
    )
  )

}
