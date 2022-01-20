find.stepfunction.chgpts = function(x){
  chg.pts = which(diff(x) != 0)
  chg.pts.plus = chg.pts+1
  keep.chg.pts = (x[chg.pts] < x[chg.pts.plus])
  c(chg.pts[keep.chg.pts], chg.pts.plus[!keep.chg.pts])
}

#' Assume that x is a histogram (or raw data?)
#' Needs to work as a standalone function
segment.fit.agnostic <- function(x,
                                 eps = 1,
                                 seed = NULL,
                                 truncated.models = FALSE,
                                 uniform.peak.threshold = 0.75,
                                 uniform.peak.stepsize = 5,
                                 remove.low.entropy = T,
                                 max.uniform = T,
                                 histogram.metric = c("jaccard", "intersection", "ks", "mse", "chisq")) {

  # Change points
  chgpts = find.stepfunction.chgpts(x)
  p.init = sort(unique(c(1, chgpts, length(x))))
  p = ftc.helen(x, p.init, eps)

  # Max Gap
  if(remove.low.entropy){
    mgaps = meaningful.gaps.local(x = x, seg.points = p, change.points = p.init)
    max.gaps = mgaps[,c("Var1", "Var2")]
    p.pairs <- remove.max.gaps.agnostic(p, max.gaps, remove.short.segment = 1)
  }

  # Fitting different models
  models = list()
  set.seed(seed)
  for(i in seq_along(p.pairs)){
    seg = p.pairs[[i]]
    seg.start = seg$start
    seg.end = seg$end
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
    best.models[['majority.vote']] <- dist.optim[[paste0("jaccard.", dist.majority)]]

    # Correcting for optimization via finding the minimum
    best.models = lapply(best.models, function(mod){
      if(mod$metric %in% c("jaccard", "intersection")){ mod$value = 1 - mod$value }
      mod
    })

    models[[i]] <- best.models
  }

  return(models)
}
