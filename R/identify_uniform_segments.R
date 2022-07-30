#' Finds the largest uniform segment that is longer than threshold
#'
#' @param x numeric vector representing the density of a histogram
#' @param metric one of `jaccard`, `intersection`, `ks`, `mse`, `chisq` indicating metrics to use for fit optimization
#' @param threshold numeric, indicating the minimum proportion of the subsegment which should be tested
#' @param stepsize integer, indicating the stepsize (relative to the histogram bins) to take in the search for the uniform subsegment
#' @param max.sd.size numeric, the number of standard deviations of the computed metric distribution away from the optimal uniform which 
#' has maximum length
#' 
#' @return A data.frame with the following columns
#' \describe{
#'     \item{a}{start index of the maximum uniform segment}
#'     \item{b}{end index of the maximum uniform segment}
#'     \item{metric}{value of the fitted metric on the segment}
#'     \item{length}{length of the segment}
#'}
#' 
#' @export
identify_uniform_segment <- function(
  x,
  metric = c("jaccard", "intersection", "ks", "mse", "chisq"),
  threshold = 0.5,
  stepsize = 1,
  max.sd.size = 1) {
  
  # Error checking
  metric <- match.arg(metric)
  
  # Set-up
  num.bins <- length(x)
  min.seg.size <- ceiling(num.bins * threshold)
  metric.func <- get(paste('histogram', metric, sep = "."))

  p.unif <- generate_uniform_distribution(x)
  res <- lapply(seq(from = 1, to = num.bins - min.seg.size, by = stepsize), function(a) {
    lapply(seq(from = min.seg.size + a, to = num.bins, by = stepsize), function(b) {
      x.sub <- x[a:b]
      p.unif.sub <- generate_uniform_distribution(x.sub)
      h.sub <- x.sub / sum(x.sub)

      m <- metric.func(h.sub, p.unif.sub)

      return(
        list('start' = a, 'end' = b, 'metric' = m)
      )
    })
  })

  res.df <- do.call(rbind.data.frame, unlist(res, recursive = F))
  res.df$length <- res.df$end - res.df$start + 1

  # Select the longest interval that is within 1 sd of the maximum
  min.metric <- min(res.df$metric, na.rm = T)
  sd.metric <- stats::sd(res.df$metric, na.rm = T)
  sd.metric <- ifelse(is.na(sd.metric), 0, sd.metric) # If there's only 1 case
  # The range in which we are looking for the minimum
  res.sd.range <- res.df[res.df$metric <= min.metric + sd.metric * max.sd.size, ]
  max.interval.index <- which.max(res.sd.range$length)
  
  return(
    res.sd.range[max.interval.index,]
  )
  
}
