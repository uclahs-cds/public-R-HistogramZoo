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

local.min = function(x) {
  which(diff(sign(diff(x)))==2)
}

local.max = function(x) {
  which(diff(sign(diff(x)))==-2)
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
  # Add end points
  m = c(1, local.min(x), length(x))
  M = local.max(x)
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

plot.segments = function(x) {
  index = seq_along(x)
  plot(x)

  tp = turnpoints(x)
  min.ind = tp$pos[tp$pits]
  max.ind = tp$pos[tp$peaks]
  points(min.ind, x[min.ind], col = "red")
  points(max.ind, x[max.ind], col = "green")
}

set.seed(13)
# Example 1, uniform
x = obs.to.int.hist(runif(10, min = 0, max = 50))
plot.segments(x)
s = ftc(x)

# Example 2, normal dist
x.norm = obs.to.int.hist(rnorm(1000) * 5)
s.norm = ftc(x.norm)
plot.segments(x.norm)
s.norm

# Example 2, messy mixture of gaussians and uniform
x.norm.mix = obs.to.int.hist(c(rnorm(100, mean = 1), rnorm(100, mean = 5), c(runif(20, min = 3, max = 5))))
s.norm.mix = ftc(x.norm.mix)
plot.segments(x.norm.mix)

