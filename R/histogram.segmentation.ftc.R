#' Take a vector of values and get the histogram for integer breaks
#' @param x the observation
#' @param add.zero.endpoints Should the left and right side be padded by 1? This will make one bin zero on each side.
#' @export
obs.to.int.hist = function(x, add.zero.endpoints = TRUE, as.df = FALSE) {
  a = floor(min(x))-1
  b = ceiling(max(x))
  breaks = a:b
  if(add.zero.endpoints) breaks = c(a - 1, breaks, b + 1)
  rtn = table(cut(x, breaks = breaks))
  names(rtn) <- breaks[2:length(breaks)]
  if(as.df) {
    rtn <- as.data.frame(rtn)
    colnames(rtn) <- c("x", "Freq")
    rtn$x <- as.integer(as.character(rtn$x))
  }
  rtn
}

# Plots the vector x of counts (or table) and the optional segment points s
plot.segments = function(x, s = NULL, threshold = 0, ...) {
  index = seq_along(x)
  if(!is.null(s)) {
    opar = par(mfrow = c(2,1), mar = c(2,2,2,2))
  }
  plot(x, type = "h", ...)

  minmax = local.minmax(x, threshold)
  min.ind = minmax$min.ind
  max.ind = minmax$max.ind
  # min.ind = find_peaks(-x, strict = FALSE)
  # max.ind = find_peaks(x, strict = TRUE)
  # min.ind = local.min(x)
  # max.ind = local.max(x)
  both.ind = intersect(min.ind, max.ind)

  # points(seq_along(x)[min.ind & !max.ind], x[min.ind & !max.ind], col = "green")
  # points(seq_along(x)[max.ind & !min.ind], x[max.ind & !min.ind], col = "red")
  # points(seq_along(x)[max.ind & min.ind], x[max.ind & min.ind], col = "orange")

  points(seq_along(x)[min.ind], x[min.ind], col = "green")
  points(seq_along(x)[max.ind], x[max.ind], col = "red")
  points(seq_along(x)[both.ind], x[both.ind], col = "orange")

  if(!is.null(s)) {
    plot(x, type = "h")
    points(s, x[s], col = "orange")
    par(opar)
  }
}

#' Kullback-Leibler divergence (Relative Entropy)
rel.entropy = function(h, p, a, b) {
  interval = a:b
  # Round to prevent floating point issues
  hab = round(sum(h[interval]), digits = 14)
  pab = round(sum(p[interval]), digits = 14)
  if(pab == 0 || pab == 1) {
    return(0)
  } else {
    hab * log(hab / pab) + (1 - hab) * log((1 - hab) / (1 - pab))
  }
}

# Estimate x with grenader estimator
grenader = function(x, increasing = T){
  if(increasing){
    est = isotone::gpava(z = 1:length(x), y = x)
    est = est$x
  } else {
    est = isotone::gpava(z = 1:length(x), y = rev(x))
    est = rev(est$x)
  }
  N = ifelse(sum(x) == 0, 1, sum(x))
  return(est / N)
}

local.min = function(x) {
  which(diff(sign(diff(x)))==2)+1
}

local.max = function(x) {
  which(diff(sign(diff(x)))==-2)+1
}

#' Get the max relative entropy in the interval
#' Computes H, the maximum H_{h,p}([a,b])
max.entropy = function(x, s = NULL, increasing = TRUE) {
  N = ifelse(sum(x) == 0, 1, sum(x))
  if(is.null(s)){s = 1:length(x)}
  L = length(s)

  # Prob distributions
  h = x/N
  p = grenader(x, increasing)

  max.rel.entropy = -Inf
  for(a in 1:L) {
    for(b in a:L) {
      max.rel.entropy = max(max.rel.entropy, rel.entropy(h, p, s[a], s[b]), na.rm = TRUE)
    }
  }
  max.rel.entropy
}

#' Finds the local minima m and maxima M such that
#' m_1 < M_1 < m_2 < M_2 < ... < M_{K - 1} < m_{k}
#' @export
local.minmax = function(x, threshold = 0) {
  x = as.numeric(x)
  stopifnot(length(x) > 1)
  # Get the first non-equal index
  n = length(x)
  min.ind = NULL
  max.ind = NULL
  #init = x[1]
  # If x = c(1,1,1,1,2,3) then j = 4, the first index before the run ends

  # Trim the left and right side
  rle_diff = rle(diff(x))
  # Return no min/max if all equal
  if(length(rle_diff$lengths) <= 1) {
    return(list(min.ind = 1, max.ind = NULL))
  }
  startIndex = 1
  endIndex = n
  if(rle_diff$values[1] == 0) {
    startIndex = rle_diff$lengths[1] + 1
  }
  if(length(rle_diff$lengths) > 2 && tail(rle_diff$values, n = 1) == 0) {
    # Remove the last equal values
    endIndex = n - tail(rle_diff$lengths, n = 1)
  }
  x.trim = x[startIndex:endIndex]
  n.trim = length(x.trim)
  # Keep track if we last appended a min/max
  min.appended = NULL
  # Check the first point for min/max
  if(x.trim[1] < x.trim[2]) {
    min.ind = 1
    min.appended = TRUE
  } else if (x.trim[1] > x.trim[2]) {
    max.ind = 1
    min.appended = FALSE
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
        min.ind = c(min.ind, i)
        min.appended = TRUE
      } else if(min.appended &&
                ((left.diff > threshold && right.diff >= 0) ||
                (left.diff >= 0 && right.diff > threshold))) {
        max.ind = c(max.ind, i)
        min.appended = FALSE
      }
    }
  }

  if(x.trim[n.trim - 1] > x.trim[n.trim]) {
    min.ind = c(min.ind, n.trim)
    min.appended = TRUE
  } else if (x.trim[n.trim - 1] < x.trim[n.trim]) {
    max.ind = c(max.ind, n.trim)
    min.appended = FALSE
  }

  list(min.ind = min.ind + (startIndex - 1), max.ind = max.ind + (startIndex - 1))
}

# Compute the monotone cost,
monotone.cost = function(x, s = NULL, eps = 1, increasing = TRUE) {
  max.rel.entropy = max.entropy(x, s, increasing)
  N = sum(x)
  L = length(x)

  max.rel.entropy * N - log(L * (L + 1) / 2 * eps)
}

#' Fine-to-Course Algorithm from Lisani & Petro 2021
#'
#' @param x a vector (or table) of counts representing the histogram
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

ftc.helen = function(x, s = NULL, eps = 1) {
  if(is.null(s)){
    # segments
    minmax = local.minmax(x)
    lmin = minmax$max.ind
    lmax = minmax$min.ind
    s = c(1, lmin, lmax, length(x) )
    s = sort(unique(s))
  }

  # Initializing
  s.fixed = s
  K = length(s)
  cost = c(-Inf)
  J = 1

  while(!all(cost > 0) & K > 2){
    # Initialize
    cost = -Inf
    # Loop through segments
    for(i in 1:(K-2)){
      # cat(K, " ", i, "\n")
      inc.int = s[i]:s[i+1]
      inc.s = s.fixed[s.fixed >= s[i] & s.fixed <= s[i+1]] - s[i] + 1
      dec.int = s[i+1]:s[i+2]
      dec.s = s.fixed[s.fixed >= s[i+1] & s.fixed <= s[i+2]] - s[i+1] + 1
      cost_i = monotone.cost(x[inc.int], eps = eps, s=inc.s, increasing = TRUE)
      cost_d = monotone.cost(x[dec.int], eps = eps, s=dec.s, increasing = FALSE)
      cost[i] = min(cost_i, cost_d)
    }
    # Removing minimum cost
    min.cost.index = which.min(cost)
    min.cost = cost[min.cost.index]
    if(length(min.cost) > 0 && min.cost < 0){
      s = s[-(min.cost.index+1)]
    }
    # Update
    K = length(s)
  }

  # Return the final list of minima
  s
}
