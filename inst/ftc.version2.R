

local.min = function(x) {
  which(diff(sign(diff(x)))==2)+1
}

local.max = function(x) {
  which(diff(sign(diff(x)))==-2)+1
}


#' Fine-to-Course Algorithm from Lisani & Petro 2021
#'
#' @param x a vector (or table) of counts representing the histogram
#' @param maxJ TODO
#' @param threshold TODO
#' @export
ftc = function(x, maxJ = Inf, threshold = 0) {
  minmax = local.minmax(x, threshold = threshold)
  # Add end points
  m = minmax$min.ind
  M = minmax$max.ind
  
  names(m) = NULL
  names(M) = NULL
  K = length(m)
  J = 1
  while(J < K && J <= maxJ) {
    costs = -Inf
    while(length(costs) > 0 && min(costs) < 0) {
      costs = NULL
      # Boolean vector indicating if the cost in increasing or decreasing
      costs.inc = NULL
      K = length(m)
      if(K - J - 1 > 0) {
        for(i in 1:(K - J - 1)) {
          inc.int = m[i]:M[i + J]
          dec.int = M[i]:m[i + J + 1]
          cost_i = monotone.cost(x[inc.int], increasing = TRUE)
          cost_d = monotone.cost(x[dec.int], increasing = FALSE)
          
          # Keep track of the minimum cost, and whether it is increasing (=TRUE) or decreasing
          costs[i] = min(cost_i, cost_d)
          costs.inc[i] = cost_i < cost_d
        }
      }
      
      min.cost.index = which.min(costs)
      if(length(min.cost.index) > 0) {
        min.cost = costs[min.cost.index]
        if(min.cost < 0) {
          # Remove the indicies
          remove.indices = seq(min.cost.index + 1, min.cost.index + J)
          m = m[-remove.indices]
          
          # If it is increasing
          if(costs.inc[min.cost.index]) {
            remove.indices.M = seq(min.cost.index, min.cost.index + J - 1)
          } else {
            # Decreasing
            remove.indices.M = seq(min.cost.index + 1, min.cost.index + J)
          }
          M = M[-remove.indices.M]
        }
      }
    }
    J = J + 1
  }
  # Return the final list of minima
  # If we selected the false start point, then increment index
  if(m[1] == 1 && x[m[1]] == 0) m[1] = 2
  L = length(m)
  n = length(x)
  # If we selected the false end point, then decrement the index
  if(m[L] == n && x[m[L]] == 0) m[L] = n - 1
  
  m
}