
#' Uniform density generation
#' @param x A vector
#'
#' @return A vector of uniform density with length of x
generate_uniform_distribution <- function(x){
  return(
    rep(1/length(x), length(x))
  )
}

#' Calculate whether h[a,b] > p[a,b]
#' @param h density of a distribution
#' @param p density of the reference distribution
#' @param a start index to subset density
#' @param b end index to subset density
#'
#' @return logical
calculate_probability_difference <- function(h, p, a, b){
  interval <- a:b
  # Round to prevent floating point issues
  hab <- round(sum(h[interval]), digits = 14)
  pab <- round(sum(p[interval]), digits = 14)
  return(hab > pab)
}

#' Determines whether h[a,b] is a meaningful interval
#' @param h density of a distribution
#' @param p density of the reference distribution
#' @param a start index to subset density
#' @param b end index to subset density
#' @param N number of data points
#' @param L length of interval
#'
#' @return a vector of length 2: binary indicating whether
#' h[a,b] is a meaningful interval, relative entropy of the interval
meaningful_interval <- function(h, p, a, b, N, L){
  rel_entropy <- relative_entropy(h, p, a, b)
  prob_diff <- calculate_probability_difference(h, p, a, b)
  return(
    c(mint = rel_entropy >= (1/N)*log(L*(L+1)/2) && prob_diff,
      entropy = rel_entropy)
  )
}

#' Determines whether h[a,b] is a meaningful gap
#' @param h density of a distribution
#' @param p density of the reference distribution
#' @param a start index to subset density
#' @param b end index to subset density
#' @param N number of data points
#' @param L length of interval
#'
#' @return a vector of length 2: binary indicating whether
#' h[a,b] is a meaningful gap, relative entropy of the interval
meaningful_gap <- function(h, p, a, b, N, L){
  rel_entropy <- relative_entropy(h, p, a, b)
  if(is.na(rel_entropy)) rel_entropy <- Inf
  prob_diff <- calculate_probability_difference(h, p, a, b)
  mgap <- (rel_entropy >= (1/N)*log(L*(L+1)/2) && !prob_diff) || (all(h == 0))
  return(
    c(mgap = mgap, entropy = rel_entropy)
  )
}

#' Identifies the maximum meaningful interval from a set of overlapping
#' intervals, can also be applied to find the max gap
#'
#' @param x a data.frame with columns
#' \describe{
#'     \item{start}{start index of the histogram count vector}
#'     \item{end}{end index of the histogram count vector}
#'     \item{entropy}{relative entropy of the segment to a uniform distribution}
#' }
#'
#' @return a data.frame with the same columns + the index of the preserved
#' rows from the original data.frame
maximal_meaningful_interval <- function(x) {
  curr.df <- x
  max.intervals <- data.frame()
  while(nrow(curr.df) > 0) {
    maximum.entropy.index <- which.max(curr.df$entropy)
    maximum.entropy <- curr.df[maximum.entropy.index, ]
    maximum.entropy.seq <- seq(maximum.entropy$start, maximum.entropy$end)
    max.intervals <- rbind(max.intervals, maximum.entropy)

    # Find all of the segments that overlap.
    # These will all be less than the maximum
    overlap.max <- mapply(function(from, to) {
      s <- seq(from, to)
      length(intersect(maximum.entropy.seq, s)) > 0
    }, from = curr.df$start, to = curr.df$end)

    curr.df = curr.df[!overlap.max, ]
  }
  # Preserve the old index
  max.intervals$index <- rownames(max.intervals)
  rownames(max.intervals) <- NULL
  return(max.intervals)
}


#' Find all meaningful gaps
#'
#' @param x histogram (vector of counts)
#' @param change.points Change points
#'
#' @return A data.frame with columns
#' \describe{
#'     \item{start}{start index of the histogram count vector}
#'     \item{end}{end index of the histogram count vector}
#'     \item{mgap}{max gap score}
#'     \item{entropy}{relative entropy of the segment to a uniform distribution}
#' }
#' @export
find_all_meaningful_gap <- function(x, change.points) {
  todo <- expand.grid(start = change.points, end = change.points)
  todo <- todo[todo$end > todo$start,]

  mgap <- do.call(rbind, lapply(1:nrow(todo), function(i) {
    meaningful_gap(
      h = x/sum(x),
      p = generate_uniform_distribution(x),
      a = todo$start[i],
      b = todo$end[i],
      N = sum(x),
      L = length(x)
    )
  }))
  df <- cbind(todo, mgap)
  # df = df[order(df$end, df$start),]
  df <- df[order(df$entropy), ]
  df$mgap <- as.numeric(df$mgap)
  # df$scaled_entropy = (df$entropy - min(df$entropy, na.rm = T)) / (max(df$entropy, na.rm = T) - min(df$entropy, na.rm = T))

  seg.gap.data <- df[df$mgap > 0 & !is.na(df$mgap), ]
  # seg.gap.data$index = 1:nrow(seg.gap.data)

  return(maximal_meaningful_interval(seg.gap.data))
}

#' Finds the meaningful gaps between the points in s
#'
#' @param x The histogram data
#' @param seg.points the segment points
#' @param change.points the change points
#' @param min.gap The minimum gap to be considered a meaningful gap
#'
#' @return A data.frame with columns
#' \describe{
#'     \item{start}{start index of the histogram count vector}
#'     \item{end}{end index of the histogram count vector}
#'     \item{mgap}{max gap score}
#'     \item{entropy}{relative entropy of the segment to a uniform distribution}
#' }
#' @export
meaningful_gaps_local <- function(x, seg.points, change.points, min.gap = 2) {

  max.gaps.list <- lapply(seq(2, length(seg.points)), function(i) {
    x.sub <- x[seg.points[i-1]:seg.points[i]]
    chg.pts <- change.points[change.points >= seg.points[i-1] & change.points <= seg.points[i]] - seg.points[i-1] + 1

    max.gaps <- find_all_meaningful_gap(x.sub, chg.pts)

    if(nrow(max.gaps) > 0) {
      max.gaps[, c('start','end')] <- max.gaps[, c('start','end')] + seg.points[i-1] - 1
      max.gaps$seg.start <- seg.points[i-1]
      max.gaps$seg.end <- seg.points[i]
      max.gaps
    }
  })

  max.gaps <- do.call(rbind.data.frame, max.gaps.list)
  # Remove gaps that are smaller than min.gap
  return(max.gaps[max.gaps$end - max.gaps$start >= min.gap, ])
}
