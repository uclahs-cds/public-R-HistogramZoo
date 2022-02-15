
#' Title
#'
#' @param coverage.model.obj TODO
#' @param histogram.count.threshold TODO
#' @param eps TODO
#' @param seed TODO
#' @param truncated.models TODO
#' @param uniform.peak.threshold TODO
#' @param uniform.peak.stepsize TODO
#' @param remove.low.entropy TODO
#' @param min.gap.size TODO
#' @param min.peak.size TODO
#' @param max.uniform TODO
#' @param histogram.metric TODO
#'
#' @return
#' @export
bulk.segment.fit = function(
  coverage.model.obj,
  histogram.count.threshold = 0,
  eps = 10^-4,
  seed = NULL,
  truncated.models = FALSE,
  uniform.peak.threshold = 0.75,
  uniform.peak.stepsize = 5,
  remove.low.entropy = T,
  min.gap.size = 2,
  min.peak.size = 1,
  max.uniform = T,
  histogram.metric = c("jaccard", "intersection", "ks", "mse", "chisq")
){

  # Coverage.model.obj
  cov = coverage.model.obj$histogram.coverage
  histogram.ids = names(cov)

  results = vector("list", length(histogram.ids))
  names(results) = histogram.ids
  # Running segmentation & distribution fitting
  for(i in histogram.ids){
    res = segment.fit.agnostic(
      x = cov[[i]],
      eps = eps,
      seed = seed,
      truncated.models = truncated.models,
      uniform.peak.threshold = uniform.peak.threshold,
      uniform.peak.stepsize = uniform.peak.stepsize,
      remove.low.entropy = remove.low.entropy,
      min.gap.size = min.gap.size,
      min.peak.size = min.peak.size,
      max.uniform = max.uniform,
      histogram.metric = histogram.metric
    )
    results[[i]] <- res
  }

  coverage.model.obj[['results']] <- results
  coverage.model.obj
}
