
#' Kullback-Leibler divergence (Relative Entropy)
#'
#' @param h density of a distribution
#' @param p density of the reference distribution
#' @param a start index to subset density
#' @param b end index to subset density
#'
#' @return numeric, relative entropy of the distributions `h` and `p` between `a` and `b`
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
#'
#' @return numeric vector, representing maximum relative entropy of all subsegments of `h`
#' (the distribution of normalizing `x`) and `p` (the Grenander estimator of x)
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
#' @return a vector of points representing the boundaries between histograms
#'
#' @examples \dontrun{
#' x = c(0, 1, 2, 3, 4, 3, 2, 1, 0, 0, 0, 1, 2, 1, 0)
#' s = c(1, 9, 11, 15)
#' ftc(x = x, s = s, eps = 1)
#' }
#' 
#' @rdname ftc
#' @export
ftc <- function(x, s = NULL, eps = 1) {
  UseMethod("ftc")
}

#' @rdname ftc
#' @exportS3Method ftc Histogram
ftc.Histogram <- function(x, s = NULL, eps = 1) {
  ftc.numeric(x$histogram_data, s = s, eps = eps)
}

#' @rdname ftc
#' @exportS3Method ftc numeric
ftc.numeric <- function(x, s = NULL, eps = 1) {

  # Error checking
  if(!is_equal_integer(s) | !all(s <= length(x) & s >= 0)){
    stop("s must be functional indices")
  }
  if(!is.numeric(eps) | !(eps > 0)){
    stop("eps must be a positive numeric")
  }

  # Generating change points if none are provided
  if(is.null(s)){
    # segments
    minmax <- find_local_optima(x)
    lmin <- minmax$max_ind
    lmax <- minmax$min_ind
    s <- c(1, lmin, lmax, length(x))
  }

  # Initializing
  s <- sort(unique(s))
  s_fixed <- s
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
      inc.s <- s_fixed[s_fixed >= s[i] & s_fixed <= s[i+1]] - s[i] + 1
      dec.int <- s[i+1]:s[i+2]
      dec.s <- s_fixed[s_fixed >= s[i+1] & s_fixed <= s[i+2]] - s[i+1] + 1
      cost_i <- monotone_cost(x[inc.int], eps = eps, s=inc.s, increasing = TRUE)
      cost_d <- monotone_cost(x[dec.int], eps = eps, s=dec.s, increasing = FALSE)
      cost[i] <- min(cost_i, cost_d)
    }
    # Removing minimum cost
    min_cost_index <- which.min(cost)
    min_cost <- cost[min_cost_index]
    if(length(min_cost) > 0 && min_cost < 0){
      s <- s[-(min_cost_index+1)]
    }
    # Update
    K <- length(s)
  }

  # Return the final list of minima
  return(s)
}
