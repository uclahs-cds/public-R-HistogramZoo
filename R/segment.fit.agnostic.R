#' Assume that x is a histogram (or raw data?)
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
    mgaps$Var1 = mgaps$Var1-1
    mgaps$Var2 = mgaps$Var2-1
    max.gaps = mgaps[,c("Var1", "Var2")]

    p <- remove.max.gaps.agnostic(p, max.gaps, remove.short.segment = 1)
  }

  # Fitting different models
  results = data.frame()
  models = list()
  set.seed(seed)
  for(i in 2:length(p)){
    # cat(i , "\n")
    # Extracting data
    seg.start = p[i - 1]
    seg.end = p[i]
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
    value.df = data.frame(
      "i" = i,
      "value" = unlist(lapply(dist.optim, `[[`, "value")),
      "dist" = unlist(lapply(dist.optim, `[[`, "dist")),
      "metric" = unlist(lapply(dist.optim, `[[`, "metric")),
      "params" = unlist(lapply(dist.optim, function(m) dput.str(m$par))),
      "seg.start" = unlist(lapply(dist.optim, `[[`, "seg.start")),
      "seg.end" = unlist(lapply(dist.optim, `[[`, "seg.end")),
      "final" = 0,
      stringsAsFactors=F
    )
    distr.vote = aggregate(value ~ metric, value.df, FUN = min)
    vote.df = merge(value.df, distr.vote, by = c("metric", "value"))
    distr.tally = table(vote.df$dist)
    best.distr = ifelse(sum(distr.tally == max(distr.tally))>1, vote.df$dist[vote.df$metric == "jaccard"], names(distr.tally)[which.max(distr.tally)])
    final.res = value.df[value.df$metric == "jaccard" & value.df$dist == best.distr,]
    final.res$final <- 1
    vote.df = rbind.data.frame(vote.df, final.res)

    mod.final = dist.optim[[paste0("jaccard.", best.distr)]]
    models[[i]] = mod.final
    results = rbind(results, vote.df)
  }

  # Correcting for optimization via finding the minimum
  results$value[results$metric %in% c("jaccard", "intersection")] <- 1 - results$value[results$metric %in% c("jaccard", "intersection")]

  return(results)
}
