#' Take a vector of values and get the histogram for integer breaks
#'
#' @param x A vector of integer observations
#' @param add.zero.endpoints Should the left and right side be padded by 1? This will make one bin zero on each side.
#' @param as.df Should a dataframe be returned? Default FALSE (e.g return a named vector)
#'
#' @return A vector representing the number observations from the minimum observed value to the maximum observed value
#'
#' @export
observations_to_histogram = function(x, add.zero.endpoints = TRUE, as.df = FALSE) {
  a = floor(min(x))-1
  b = ceiling(max(x))
  breaks = a:b
  if(add.zero.endpoints) breaks = c(a - 1, breaks, b + 1)
  rtn = table(cut(x, breaks = breaks))
  names(rtn) <- breaks[2:length(breaks)]
  # Convert from table to named vector
  rtn <- unlist(as.list(rtn))
  if(as.df) {
    rtn <- as.data.frame(rtn)
    colnames(rtn) <- c("x", "Freq")
    rtn$x <- as.integer(as.character(rtn$x))
  }
  rtn
}

#' Kullback-Leibler divergence (Relative Entropy)
#'
#' @param h TODO
#' @param p TODO
#' @param a interval left endpoint
#' @param b interval right endpoint
relative_entropy = function(h, p, a, b) {
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

#' Get the max relative entropy in the interval
#' Computes H, the maximum H_{h,p}([a,b])
#'
#' @param x numeric vector of counts representing a histogram
#' @param s TODO
#' @param increasing Should the Grenader estimator be increasing or descreasing?
#'
#' @export
maximum_entropy = function(x, s = NULL, increasing = TRUE) {
  N = ifelse(sum(x) == 0, 1, sum(x))
  if(is.null(s)){s = 1:length(x)}
  L = length(s)

  # Prob distributions
  h = x/N
  p = grenader(x, increasing)

  max.relative_entropy = -Inf
  for(a in 1:L) {
    for(b in a:L) {
      max.relative_entropy = max(max.relative_entropy, relative_entropy(h, p, s[a], s[b]), na.rm = TRUE)
    }
  }
  max.relative_entropy
}

# Compute the monotone cost,
monotone_cost = function(x, s = NULL, eps = 1, increasing = TRUE) {
  max.relative_entropy = maximum_entropy(x, s, increasing)
  N = sum(x)
  L = length(x)

  max.relative_entropy * N - log(L * (L + 1) / 2 * eps)
}

#' Fine-to-coarse segmentation algorithm
#'
#' @param x A vector of counts representing the bins of a histogram
#' @param s Initialization points
#' @param eps The hyperparameter epsilon which can be used to tune the scale of segmentation
#'
#' @return A vector of points representing the boundaries between histograms
#'
#' @examples \dontrun{
#' x = c(0, 1, 2, 3, 4, 3, 2, 1, 0, 0, 0, 1, 2, 1, 0)
#' s = c(1, 9, 11, 15)
#' ftc(x = x, s = s, eps = 1)
#' }
#'
#' @export
ftc = function(x, s = NULL, eps = 1) {
  if(is.null(s)){
    # segments
    minmax = local.minmax(x)
    lmin = minmax$max.ind
    lmax = minmax$min.ind
    s = c(1, lmin, lmax, length(x) )
    s = sort(unique(s))
  }

  # Initializing
  s = sort(unique(s))
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
      cost_i = monotone_cost(x[inc.int], eps = eps, s=inc.s, increasing = TRUE)
      cost_d = monotone_cost(x[dec.int], eps = eps, s=dec.s, increasing = FALSE)
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
  return(s)
}
