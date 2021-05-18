# Take a vector of values and get the histogram for integer breaks
obs.to.int.hist = function(x) {
  table(cut(x, breaks = floor(min(x)):ceiling(max(x))))
}

rel.entropy = function(h, p, a, b) {
  interval = a:b
  hab = sum(h[interval])
  pab = sum(p[interval])
  if(hab == pab){
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

# Get the max relative entropy in the interval
max.entropy = function(x, increasing = TRUE) {
  N = ifelse(sum(x) == 0, 1, sum(x))
  L = length(x)

  # Prob distributions
  h = x/N
  p = grenader(x, increasing)

  max.rel.entropy = -Inf
  for(a in 1:L) {
    for(b in a:L) {
      max.rel.entropy = max(max.rel.entropy, rel.entropy(h, p, a, b), na.rm = TRUE)
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

  # segments
  lmin = which(find_peaks(-x, strict = F))
  lmax = which(find_peaks(x, strict = F))
  s = c(1, lmin, lmax, length(x) )
  s = sort(unique(s))

  # Initializing
  K = length(s)
  cost = c(-Inf)

  while(!all(cost > 0) & K > 2){
    # Initialize
    cost = c(-Inf)
    # Loop through segments
    for(i in 1:(K-2)){
      # cat(i, " ", K , "\n")
      inc.int = s[i]:s[i+1]
      dec.int = s[i+1]:s[i+2]
      cost_i = monotone.cost(x[inc.int], increasing = TRUE)
      cost_d = monotone.cost(x[dec.int], increasing = FALSE)
      cost[i] = min(cost_i, cost_d)
    }
    # Removing minimum cost
    min.cost.index = which.min(cost)
    min.cost = cost[min.cost.index]
    if(min.cost < 0){
      s = s[-(min.cost.index+1)]
    }
    # Update
    K = length(s)
  }

  # Return the final list of minima
  s
}
