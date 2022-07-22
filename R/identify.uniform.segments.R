#' Finds the largest uniform segment that is longer than threshold
#'
#' @param x TODO
#' @param metric TODO
#' @param threshold TODO
#' @param step.size TODO
#' @param max.sd.size TODO
#' @export
identify_uniform_segment = function(
  x,
  metric = c("jaccard", "intersection", "ks", "mse", "chisq"),
  threshold = 0.5,
  step.size = 1,
  max.sd.size = 1) {
  metric = match.arg(metric)
  num.bins = length(x)
  min.seg.size = ceiling(num.bins * threshold)
  metric.func = get(paste('histogram', metric, sep = "."))

  p.unif = generate_uniform_distribution(x)
  res = lapply(seq(from = 1, to = num.bins - min.seg.size, by = step.size), function(a) {
    lapply(seq(from = min.seg.size + a, to = num.bins, by = step.size), function(b) {
      x.sub = x[a:b]
      p.unif.sub = generate_uniform_distribution(x.sub)
      h.sub = x.sub / sum(x.sub)

      m = metric.func(h.sub, p.unif.sub)

      list(a = a, b = b, metric = m)
    })
  })

  res.df = do.call(rbind.data.frame, unlist(res, recursive = F))
  colnames(res.df) = c("a", "b", "metric")
  res.df$length = res.df$b - res.df$a

  # Select the longest interval that is within 1 sd of the maximum
  min.metric = min(res.df$metric, na.rm = T)
  sd.metric = stats::sd(res.df$metric, na.rm = T)
  sd.metric = ifelse(is.na(sd.metric), 0, sd.metric) # If there's only 1 case
  # The range in which we are looking for the minimum
  res.sd.range = res.df[res.df$metric <= min.metric + sd.metric * max.sd.size, ]
  min.interval.index = which.min(res.sd.range$length)

  as.list(res.sd.range[min.interval.index,])
}
