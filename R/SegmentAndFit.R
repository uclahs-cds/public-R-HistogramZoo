find.stepfunction.chgpts = function(x){
  chg.pts = which(diff(x) != 0)
  chg.pts.plus = chg.pts+1
  keep.chg.pts = (x[chg.pts] < x[chg.pts.plus])
  c(chg.pts[keep.chg.pts], chg.pts.plus[!keep.chg.pts])
}

#' Returns the indices for consecutive elements of a vector that are greater than a specified threshold
#'
#' @param x A numeric vector
#' @param threshold numeric threshold
#'
#' @return TODO
#' @export
#'
#' @examples
#' find.consecutive.threshold(c(0,0,0,1,1,1,0,0,0,1,1,1,0,0))
#' find.consecutive.threshold(c(0,0,1,2,2,0,1,1,1,0,0), threshold = 1)
find.consecutive.threshold = function(
  x,
  threshold = 0){
  x.thresh = rle(x > threshold)
  end.coords = cumsum(x.thresh$lengths)
  start.coords = end.coords - x.thresh$lengths + 1
  start.coords.above.threshold = start.coords[x.thresh$values]
  end.coords.above.threshold = end.coords[x.thresh$values]

  list(start = start.coords.above.threshold, end = end.coords.above.threshold)
}

#' SegmentAndFit
#'
#' @param x a Histogram or HistogramList object
#' @param histogram.count.threshold TODO
#' @param eps TODO
#' @param seed TODO
#' @param truncated.models TODO
#' @param uniform.peak.threshold TODO
#' @param uniform.peak.stepsize TODO
#' @param remove.low.entropy TODO
#' @param min.gap.size TODO
#' @param max.uniform TODO
#' @param histogram.metric TODO
#' @param min.peak.size TODO
#'
#' @return TODO
#' @export
#' @example \dontrun{
#' x = Histogram(c(0, 0, 1, 2, 3, 2, 1, 2, 3, 4, 5, 3, 1, 0))
#' res = SegmentAndFit(x)
#' }
SegmentAndFit <- function(
  x,
  histogram.count.threshold = 0,
  eps = 1,
  seed = NULL,
  truncated.models = FALSE,
  uniform.peak.threshold = 0.75,
  uniform.peak.stepsize = 5,
  remove.low.entropy = T,
  min.gap.size = 2,
  min.peak.size = 2,
  max.uniform = T,
  histogram.metric = c("jaccard", "intersection", "ks", "mse", "chisq")
  ) {

  # Checking types
  stopifnot(inherits(x, "Histogram"))
  histogram.metric = match.arg(histogram.metric, several.ok = T)

  # Change points
  chgpts = find.stepfunction.chgpts(x)
  # Looking for regions that surpass a hard count threshold
  x.segs = as.data.frame(find.consecutive.threshold(x, threshold = histogram.count.threshold))
  x.segs = x.segs[x.segs$start != x.segs$end,]

  all.points = apply(x.segs, 1, function(segs) {
    p.init = unname(c(segs['start'], chgpts[chgpts > segs['start'] & chgpts < segs['end']], segs['end']))
    p.init = sort(unique(p.init)) # meaningful gaps local also needs p.init to be sorted so temporarily adding this back
    p = ftc(x, p.init, eps)

    # Max Gap
    if(remove.low.entropy) {
      mgaps =  meaningful.gaps.local(x = x, seg.points = p, change.points = p.init, min.gap = min.gap.size)
      p = p[(abs(p - segs['start']) > min.peak.size & abs(p - segs['end']) > min.peak.size) | p %in% segs]
      p.pairs = remove.max.gaps.agnostic(p = p, max.gaps = mgaps, remove.short.segment = min.peak.size) # remove.short.segment can also be used to filter min.peak.size, but doesn't extend to non remove low entropy cases
    } else {
      p = p[(abs(p - segs['start']) > min.gap.size & abs(p - segs['end']) > min.peak.size) | p %in% segs]
      p.pairs = index.to.start.end(p)
    }

    p.pairs
  })

  # Combine the results from each segment
  all.points = do.call('rbind.data.frame', all.points)
  all.points = all.points[(all.points$end - all.points$start + 1) > min.peak.size,, drop = FALSE] # Current solution for min.peak.size, open to alternatives
  rownames(all.points) <- NULL # use reset.rownames?

  # Fitting different models
  models = list()
  set.seed(seed)
  for(i in seq_len(nrow(all.points))) {
    seg = all.points[i, ]
    seg.start = seg[['start']]
    seg.end = seg[['end']]
    seg.len = seg.end - seg.start + 1
    bin.data = x[seg.start:seg.end]

    dist.optim = fit.distributions.optim(bin.data, metric = histogram.metric, truncated = truncated.models)
    dist.optim = lapply(dist.optim, function(y) {
      y$seg.start = seg.start
      y$seg.end = seg.end
      y
    })

    # Find the maximum uniform segment
    if(max.uniform & seg.len > uniform.peak.stepsize & seg.len > ceiling(uniform.peak.threshold*seg.len)){
      max.unif.results = find.uniform.segment(bin.data, metric = histogram.metric, threshold = uniform.peak.threshold, step.size = uniform.peak.stepsize, max.sd.size = 0)
      # Use the maximum segment
      unif.segment = unlist(max.unif.results[c('a', 'b')])
      bin.data.subset = bin.data[unif.segment[1]:unif.segment[2]]
      # Fit uniform distribution on maximum uniform segment
      dist.optim.subset = fit.distributions.optim(bin.data.subset, metric = histogram.metric, truncated = FALSE, distr = "unif")
      # Adjust the segment starts from the shifted max uniform segment
      dist.optim.subset = lapply(dist.optim.subset, function(y) {
        y$seg.start = unif.segment[1] + seg.start - 1
        y$seg.end = unif.segment[2] + seg.start - 1
        y$dist = "unif"
        y
      })
      for(munif in names(dist.optim.subset)){dist.optim[[munif]] <- dist.optim.subset[[munif]]}
    }

    # Metric Voting
    metric.idx = sapply(dist.optim, `[[`, "metric")
    best.models = lapply(histogram.metric, function(met){
      idx = which(metric.idx == met)
      dist.optim.met = dist.optim[idx]
      val = lapply(dist.optim.met, `[[`, "value")
      dist.optim.met[which.min(val)]
    })
    best.models = unlist(best.models, recursive = F)
    dist.extract = sapply(best.models, `[[`, "dist")
    dist.vote = sort(table(dist.extract), decreasing = T)
    if(length(dist.vote) > 1 & dist.vote[1] == dist.vote[2]){
       mets = sapply(best.models, `[[`, "metric")
       dist.majority = dist.extract[which(mets == "jaccard")]
    } else {
      dist.majority = names(dist.vote)[1]
    }
    best.models[['majority.vote']] = dist.optim[[paste0("jaccard.", dist.majority)]]

    # Correcting for optimization via finding the minimum
    best.models = lapply(best.models, function(mod){
      if(mod$metric %in% c("jaccard", "intersection")){ mod$value = 1 - mod$value }
      mod
    })

    best.models[['majority.vote']]$metric = "Consensus"

    models[[i]] <- best.models
  }

  attr(x, "models") <- models
  attr(x, "p") <- all.points
  # class(x) <- c("HistogramFit", class(x))

  return(x)
}
