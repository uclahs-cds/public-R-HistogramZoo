
#' Kullback-Leibler divergence (Relative Entropy)
#'
#' @param h density of a distribution
#' @param p density of the reference distribution
#' @param a start index to subset density
#' @param b end index to subset density
relative_entropy <- function(h, p, a, b) {
  interval <- a:b
  # Round to prevent floating point issues
  hab <- round(sum(h[interval]), digits = 14)
  pab <- round(sum(p[interval]), digits = 14)
  if(pab == 0 || pab == 1) {
    return(0)
  } else {
    return(
      hab * log(hab / pab) + (1 - hab) * log((1 - hab) / (1 - pab))
    )
  }
}

#' Calculate the Grenander estimator of a given density vector
#'
#' @param x numeric vector representing the density of a histogram
#' @param increasing logical, whether the Grenander estimator should assume 
#' x is increasing or decreasing 
#'
#' @return numeric vector with the same length as `x` representing the 
#' Grenander estimator of `x`
grenader <- function(x, increasing = T){
  if(increasing){
    est <- isotone::gpava(z = 1:length(x), y = x)
    est <- est$x
  } else {
    est <- isotone::gpava(z = 1:length(x), y = rev(x))
    est <- rev(est$x)
  }
  N <- ifelse(sum(x) == 0, 1, sum(x))
  return(est / N)
}

#' Computes H, the maximum H_{h,p}([a,b])
#'
#' @param x numeric vector representing the density of a histogram
#' @param s  numeric vector representing initial segmentation indices
#' @param increasing logical, whether the Grenander estimator should assume 
#' x is increasing or decreasing 
maximum_entropy <- function(x, s = NULL, increasing = TRUE) {
  N <- ifelse(sum(x) == 0, 1, sum(x))
  if(is.null(s)){s <- 1:length(x)}
  L <- length(s)

  # Prob distributions
  h <- x/N
  p <- grenader(x, increasing)

  max_relative_entropy = -Inf
  for(a in 1:L) {
    for(b in a:L) {
      max_relative_entropy <- max(max_relative_entropy, relative_entropy(h, p, s[a], s[b]), na.rm = TRUE)
    }
  }
  return(
    max_relative_entropy
  )
}

#' Compute the monotone cost,
#'
#' @param x numeric vector representing the density of a histogram
#' @param s numeric vector representing initial segmentation indices
#' @param eps numeric hyperparameter epsilon used to scale the resolution of segmentation
#' @param increasing logical, whether the Grenander estimator should assume 
#' x is increasing or decreasing 
#' 
#' @return numeric, monotone cost
monotone_cost <- function(x, s = NULL, eps = 1, increasing = TRUE) {
  max_relative_entropy <- maximum_entropy(x, s, increasing)
  N <- sum(x)
  L <- length(x)

  return(
    max_relative_entropy * N - log(L * (L + 1) / 2 * eps)
  )
}

#' Fine-to-coarse segmentation algorithm
#'
#' @param x numeric vector representing the density of a histogram
#' @param s numeric vector representing initial segmentation indices
#' @param eps numeric hyperparameter epsilon used to scale the resolution of segmentation 
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
ftc <- function(x, s = NULL, eps = 1) {
  
  # Error checking
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(s)) # Technically, we want to check that s is an index
  stopifnot(all(s <= length(x) & s >= 0))
  stopifnot(is.numeric(eps)) # Other criteria on eps?
  
  # Generating change points if none are provided
  if(is.null(s)){
    # segments
    minmax <- local_min_max(x)
    lmin <- minmax$max.ind
    lmax <- minmax$min.ind
    s <- c(1, lmin, lmax, length(x) )
    s <- sort(unique(s))
  }
  
  # Initializing
  s <- sort(unique(s))
  s.fixed <- s
  K <- length(s)
  cost <- c(-Inf)
  J <- 1

  while(!all(cost > 0) & K > 2){
    # Initialize
    cost <- -Inf
    # Loop through segments
    for(i in 1:(K-2)){
      # cat(K, " ", i, "\n")
      inc.int <- s[i]:s[i+1]
      inc.s <- s.fixed[s.fixed >= s[i] & s.fixed <= s[i+1]] - s[i] + 1
      dec.int <- s[i+1]:s[i+2]
      dec.s <- s.fixed[s.fixed >= s[i+1] & s.fixed <= s[i+2]] - s[i+1] + 1
      cost_i <- monotone_cost(x[inc.int], eps = eps, s=inc.s, increasing = TRUE)
      cost_d <- monotone_cost(x[dec.int], eps = eps, s=dec.s, increasing = FALSE)
      cost[i] <- min(cost_i, cost_d)
    }
    # Removing minimum cost
    min.cost.index <- which.min(cost)
    min.cost <- cost[min.cost.index]
    if(length(min.cost) > 0 && min.cost < 0){
      s <- s[-(min.cost.index+1)]
    }
    # Update
    K <- length(s)
  }

  # Return the final list of minima
  return(s)
}
