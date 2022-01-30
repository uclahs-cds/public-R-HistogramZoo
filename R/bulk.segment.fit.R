
#' Title
#'
#' @param coverage.model.obj
#' @param eps
#' @param seed
#' @param truncated.models
#' @param uniform.peak.threshold
#' @param uniform.peak.stepsize
#' @param remove.low.entropy
#' @param min.peak.size TODO
#' @param max.uniform
#' @param histogram.metric
#'
#' @return
#' @export
#'
#' @examples
bulk.segment.fit = function(
  coverage.model.obj,
  histogram.count.threshold = 1,
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

#' Takes a GRanges object and creates a BED12 GRanges object
#'
#' @param gr
#'
#' @return
#' @export
#'
#' @examples
#' TODO: Check that this is compatible with Base 0 and Base 1 formatting
bed6tobed12 = function(
  gr
){
  bed12.gr = range(gr)
  # blockCount
  S4Vectors::mcols(bed12.gr)$blockCount = length(gr)
  # blockSizes
  blockSizes = IRanges::width(gr)
  S4Vectors::mcols(bed12.gr)$blockSizes = paste0(blockSizes, ",", collapse = "")
  # blockStarts
  blockStarts = GenomicRanges::start(gr) - GenomicRanges::start(bed12.gr)
  S4Vectors::mcols(bed12.gr)$blockStarts = paste0(blockStarts, ",", collapse = "")
  # Return GRanges object
  bed12.gr
}

#' Formats results
#'
#' @param coverage.model.obj
#' @param output.format
#'
#' @return
#' @export
#'
#' @examples
#' TODO: Create a reasonable format for basic histograms
summarize.results = function(
  coverage.model.obj,
  output.format = c("stats.only", "BED6", "BED12")
){

  output.format = match.arg(output.format)

  # Coverage.model.obj
  gene.model = coverage.model.obj$gene.model
  histogram.ids = names(gene.model)
  histogram.bin.size = coverage.model.obj$histogram.bin.size
  results = lapply(coverage.model.obj$results, `[[`, "models")

  formatted.results = GenomicRanges::GRanges()
  for(i in histogram.ids){
    x = gene.model[i]
    x = unlist(x)
    bins = GenomicRanges::tile(x = x, width = histogram.bin.size)
    bins = unlist(bins)
    res = results[[i]]
    for(j in seq_along(res)){
      id =  paste0(i, ":", j)
      final.mod = res[[j]]$majority.vote
      f.res = GenomicRanges::reduce(bins[final.mod$seg.start:final.mod$seg.end])
      f.res = if(output.format == "BED12") bed6tobed12(f.res)
      S4Vectors::mcols(f.res)$name = id
      S4Vectors::mcols(f.res)[,c( "value", "metric", "dist")] = final.mod[c( "value", "metric", "dist")]
      S4Vectors::mcols(f.res)$params = dput.str(final.mod$par)
      formatted.results[id] = f.res
    }
  }
  # Reordering the columns
  S4Vectors::mcols(formatted.results) = S4Vectors::mcols(formatted.results)[,c("name", "value", "metric", "dist", "params", "blockCount", "blockSizes", "blockStarts")]
  data.frame(formatted.results)
}
