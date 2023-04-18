
#' Uniform density generation
#'
#' @param x a numeric vector representing the density of a histogram
#'
#' @return a vector of uniform density with length of x
generate_uniform_distribution <- function(x){
  return(
    rep(1/length(x), length(x))
  )
}

#' Calculate whether h[a,b] > p[a,b]
#'
#' @param h density of a distribution
#' @param p density of the reference distribution
#' @param a start index to subset density
#' @param b end index to subset density
#'
#' @return logical, if h[a,b] > p[a,b]
calculate_probability_difference <- function(h, p, a, b){
  interval <- a:b
  # Round to prevent floating point issues
  hab <- round(sum(h[interval]), digits = 14)
  pab <- round(sum(p[interval]), digits = 14)
  return(hab > pab)
}

#' Determines whether h[a,b] is a meaningful interval
#'
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
#'
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
  curr_df <- x
  max_intervals <- data.frame()
  while(nrow(curr_df) > 0) {
    maximum_entropy_index <- which.max(curr_df$entropy)
    maximum_entropy <- curr_df[maximum_entropy_index, ]
    maximum_entropy_seq <- seq(maximum_entropy$start, maximum_entropy$end)
    max_intervals <- rbind(max_intervals, maximum_entropy)

    # Find all of the segments that overlap.
    # These will all be less than the maximum
    overlap_max <- mapply(function(from, to) {
      s <- seq(from, to)
      length(intersect(maximum_entropy_seq, s)) > 0
    }, from = curr_df$start, to = curr_df$end)

    curr_df = curr_df[!overlap_max, ]
  }
  # Preserve the old index
  max_intervals$index <- rownames(max_intervals)
  rownames(max_intervals) <- NULL
  return(max_intervals)
}


#' Find all meaningful gaps
#'
#' @param x a numeric vector representing the density of a histogram
#' @param change_points the change points (e.g. local min/max) in the vector x
#'
#' @return A data.frame with columns
#' \describe{
#'     \item{start}{start index of the histogram count vector}
#'     \item{end}{end index of the histogram count vector}
#'     \item{mgap}{max gap score}
#'     \item{entropy}{relative entropy of the segment to a uniform distribution}
#' }
#' @export
find_all_meaningful_gap <- function(x, change_points) {
  UseMethod('find_all_meaningful_gap')
}

#' @exportS3Method find_all_meaningful_gap Histogram
find_all_meaningful_gap.Histogram <- function(x, change_points){
  find_all_meaningful_gap.numeric(x = x$histogram_data, change_points = change_points)
}

#' @exportS3Method find_all_meaningful_gap numeric
find_all_meaningful_gap.numeric <- function(x, change_points){
  
  # Error checking
  if(!is_equal_integer(change_points) | !all(change_points <= length(x) & change_points >= 0)){
    stop("change_points must be functional indices")
  }

  # Generating a todo set of segments
  todo <- expand.grid(start = change_points, end = change_points)
  todo <- todo[todo$end > todo$start,]

  # Identifying meaningful gaps
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

  seg_gap_data <- df[df$mgap > 0 & !is.na(df$mgap), ]
  # seg_gap_data$index = 1:nrow(seg_gap_data)

  return(maximal_meaningful_interval(seg_gap_data))
}

#' Finds the meaningful gaps between the points in s
#'
#' @param x a numeric vector representing the density of a histogram
#' @param seg_points the segment points in the vector x
#' @param change_points the change points (e.g. local min/max) in the vector x
#' @param min_gap The minimum gap to be considered a meaningful gap
#'
#' @return a data.frame with columns
#' \describe{
#'     \item{start}{start index of the histogram count vector}
#'     \item{end}{end index of the histogram count vector}
#'     \item{mgap}{max gap score}
#'     \item{entropy}{relative entropy of the segment to a uniform distribution}
#' }
#' @export
meaningful_gaps_local <- function(x, seg_points, change_points, min_gap = 2) {
  UseMethod('meaningful_gaps_local')
}

#' @exportS3Method meaningful_gaps_local Histogram
meaningful_gaps_local.Histogram <- function(x, seg_points, change_points, min_gap = 2) {
  meaningful_gaps_local.numeric(x = x$histogram_data, seg_points = seg_points, change_points = change_points, min_gap = min_gap)
}

#' @exportS3Method meaningful_gaps_local numeric
meaningful_gaps_local.numeric <- function(x, seg_points, change_points, min_gap = 2) {

  # Error checking
  if(!is_equal_integer(seg_points) | !all(seg_points <= length(x) & seg_points >= 0)){
    stop("seg_points must be functional indices")
  }
  if(!is_equal_integer(change_points) | !all(change_points <= length(x) & change_points >= 0)){
    stop("change_points must be functional indices")
  }
  if(!is_equal_integer(min_gap) | length(min_gap) != 1){
    stop("min_gap must be a numeric integer of length 1")
  }

  seg_points <- sort(unique(seg_points))
  change_points <- sort(unique(change_points))

  # Identifying max gaps
  max_gaps_list <- lapply(seq(2, length(seg_points)), function(i) {

    x_sub <- x[seg_points[i-1]:seg_points[i]]
    chg_pts <- change_points[change_points >= seg_points[i-1] & change_points <= seg_points[i]] - seg_points[i-1] + 1
    chg_pts <- unique(c(1, chg_pts, length(x_sub)))

    max_gaps <- find_all_meaningful_gap(x_sub, chg_pts)

    if(nrow(max_gaps) > 0) {
      max_gaps[, c('start','end')] <- max_gaps[, c('start','end')] + seg_points[i-1] - 1
      max_gaps$seg_start <- seg_points[i-1]
      max_gaps$seg_end <- seg_points[i]
      max_gaps
    }
  })

  max_gaps <- do.call(rbind.data.frame, max_gaps_list)
  # Remove gaps that are smaller than min_gap
  return(max_gaps[max_gaps$end - max_gaps$start >= min_gap, ])
}
