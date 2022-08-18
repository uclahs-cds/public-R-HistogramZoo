

find_local_optima_draft <- function(x, threshold = 0, flat_endpoints = T){
  
  # Error checking
  x <- as.numeric(x)
  stopifnot(length(x) > 1)
  
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
  x_sign <- sign(x_change)
  x_min <- which(diff(x_sign) == 2)+1
  x_max <- which(diff(x_sign) == -2)+1
  
  # If flat endpoints, take both ends of the flat segments as indices
  if(flat_endpoints){
    x_min <- c(unflat_idx[x_min], unflat_idx[x_min - 1] + 1)
    x_max <- c(unflat_idx[x_max], unflat_idx[x_max - 1] + 1)
  } else { # Take the midpoint
    x_min <- rowMeans(cbind(unflat_idx[x_min], unflat_idx[x_min - 1] + 1))
    x_min <- round(x_min)
    x_max <- rowMeans(cbind(unflat_idx[x_max], unflat_idx[x_max - 1] + 1))
    x_max <- round(x_max)
  }
  
  # Start & end points
  start_point <- unique(c(1, unflat_idx[1]))
  end_point <- unique(c(length(x), unflat_idx[length_unflat-1]+1))
  if(!flat_endpoints) {
    start_point <- round(mean(start_point))
    end_point <- round(mean(end_point))
  }
  
  # Assign start & end points to min/max
  if(unflat_x$values[2] > unflat_x$values[1]){
    x_min <- c(start_point, x_min)
  } else {
    x_max <- c(start_point, x_max)
  }
  
  if(unflat_x$values[length_unflat] < unflat_x$values[length_unflat - 1]){
    x_min <- c(x_min, end_point)
  } else {
    x_max <- c(x_max, end_point)
  }
  
  # Final sort
  x_min <- sort(unique(x_min))
  x_max <- sort(unique(x_max))
  
  # Return results
  return( 
    list('min_ind' = x_min, 'max_ind' = x_max)
  )
  
}