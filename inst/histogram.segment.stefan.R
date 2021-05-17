# Take a vector of values and get the histogram for integer breaks
obs.to.int.hist = function(x) {
  table(cut(x, breaks = floor(min(x)):ceiling(max(x))))
}

rel.entropy = function(h, p, a, b) {
  interval = a:b
  hab = sum(h[interval]) / sum(h)
  pab = sum(p[interval])
  hab * log(hab / pab) + (1 - hab) * log((1 - hab) / (1 - pab))
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

  return(est / sum(est))
}

local.minmax = function(x) {
  n <- length(x)
  min.bool = find_peaks(-x, strict = FALSE)
  max.bool = find_peaks(x, strict = FALSE)
  rle_valleys <- rle(as.logical(min.bool))
  rle_peaks <- rle(as.logical(max.bool))

  minima_lengths <- rle_valleys$lengths[rle_valleys$values]
  minima_end_indices <- cumsum(rle_valleys$lengths)[rle_valleys$values]
  minima_start_indices <- minima_end_indices - minima_lengths + 1

  maxima_lengths <- rle_peaks$lengths[rle_peaks$values]
  maxima_end_indices <- cumsum(rle_peaks$lengths)[rle_peaks$values]
  maxima_start_indices <- maxima_end_indices - maxima_lengths + 1

  i.min <- minima_start_indices[minima_lengths == 2]
  i.max <- maxima_start_indices[maxima_lengths == 2]

  # Remove the local minima after the duplicated
  min.bool[i.min + 1] <- FALSE
  max.bool[i.max + 1] <- FALSE

  valley_runs = minima_lengths > 2
  valley_lengths = minima_lengths[valley_runs]
  valley_start = minima_start_indices[valley_runs]
  for(j in seq_len(sum(valley_runs))) {
    start.i <- valley_start[j]
    run_length <- valley_lengths[j]

    # Remove all the internal local maxima
    # Keep the end points since we need these to be local maxima, since
    #   the next points will be local maxima
    min.bool[(start.i + 1):(start.i + run_length - 2)] <- FALSE
    max.bool[(start.i + 1):(start.i + run_length - 2)] <- FALSE
    # Add a single local maxima point in between
    max.i <- floor((2 * start.i + run_length - 1) / 2)
    max.bool[max.i] <- TRUE
    stopifnot(start.i < max.i && max.i < start.i + run_length - 1)
  }

  peak_runs = maxima_lengths > 2
  peak_lengths = maxima_lengths[peak_runs]
  peak_start = maxima_start_indices[peak_runs]
  for(j in seq_len(sum(peak_runs))) {
    start.i <- peak_start[j]
    run_length <- peak_lengths[j]

    # Remove all the internal local minima
    # Keep the end points since we need these to be local minima, since
    #   the next points will be local maxima
    min.bool[(start.i + 1):(start.i + run_length - 2)] <- FALSE
    max.bool[(start.i + 1):(start.i + run_length - 2)] <- FALSE
    # Add a single local maxima point in between
    min.i <- floor((2 * start.i + run_length - 1) / 2)
    min.bool[min.i] <- TRUE
    stopifnot(start.i < min.i && min.i < start.i + run_length - 1)
  }

  # Check that we don't have a min/max as the same point
  if(max(as.integer(min.bool) + as.integer(max.bool)) > 1) {
    browser()
    stop("Overlapping min/max at the same point")
  }
  # stopifnot()
  # Check that we alternate min and max
  check.vec <- rep("none", length(min.bool))
  check.vec[min.bool] <- "min"
  check.vec[max.bool] <- "max"
  check.rle <- rle(check.vec)
  check.minmax.lengths <- check.rle$lengths[check.rle$values %in% c("min", "max")]
  stopifnot(max(check.minmax.lengths) <= 1)

  min.ind = seq_along(x)[min.bool]
  max.ind = seq_along(x)[max.bool]

  # Start with local maximum, add first point to minima
  if(max.ind[1] < min.ind[1]) {
    min.ind <- c(1, min.ind)
  } else {
    max.ind <- c(1, max.ind)
  }

  n.max <- length(max.ind)
  n.min <- length(min.ind)
  # End with local maximum, add last point as local minimum
  if(max.ind[n.max] > min.ind[n.min]) {
    min.ind <- c(min.ind, n)
  } else {
    max.ind <- c(max.ind, n)
  }

  list(min.ind = min.ind, max.ind = max.ind)
}

# Get the max relative entropy in the interval
max.entropy = function(x, increasing = TRUE) {
  N = sum(x)
  L = length(x)
  p = grenader(x, increasing)

  max.rel.entropy = -Inf
  for(a in 1:L) {
    for(b in a:L) {
      max.rel.entropy = max(max.rel.entropy, rel.entropy(x, p, a, b), na.rm = TRUE)
    }
  }
  max.rel.entropy
}

monotone.cost = function(x, eps = 1, increasing = TRUE) {
  max.rel.entropy = max.entropy(x, increasing)
  N = sum(x)
  L = length(x)

  max.rel.entropy * N - log(L * (L + 1) / 2 * eps)
}

ftc = function(x) {
  minmax = local.minmax(x)
  # Add end points
  #m = c(1, local.min(x), length(x))
  #M = local.max(x)
  m = minmax$min.ind
  M = minmax$max.ind

  if(m[1] > M[1]) {
    minmax = local.minmax(-x)
    m = minmax$min.ind
    M = minmax$max.ind
  }

  names(m) = NULL
  names(M) = NULL
  K = length(m)
  J = 1
  while(J < K) {
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
  m
}

plot.segments = function(x, s = NULL) {
  index = seq_along(x)
  if(!is.null(s)) {
    opar = par(mfrow = c(2,1), mar = c(2,2,2,2))
  }
  plot(x)

  minmax = local.minmax(x)
  min.ind = minmax$min.ind
  max.ind = minmax$max.ind

  points(seq_along(x)[min.ind], x[min.ind], col = "green")
  points(seq_along(x)[max.ind], x[max.ind], col = "red")

  if(!is.null(s)) {
    plot(x)
    points(s, x[s], col = "orange")
    par(opar)
  }
}

set.seed(13)
# Example 1, uniform
x = obs.to.int.hist(runif(1000, min = 0, max = 50))
s = ftc(x)
plot.segments(x, s)

#set.seed(13)
# set.seed(26)
# Example 2, normal dist
x.norm = obs.to.int.hist(rnorm(1000) * 5)
s.norm = ftc(x.norm)
plot.segments(x.norm, s.norm)

# Example 2, messy mixture of gaussians and uniform
x.norm.mix = obs.to.int.hist(c(rnorm(25, mean = 1), rnorm(10, mean = 5), c(runif(20, min = 3, max = 5))))
s.norm.mix = ftc(x.norm.mix)
plot.segments(x.norm.mix, s.norm.mix)

