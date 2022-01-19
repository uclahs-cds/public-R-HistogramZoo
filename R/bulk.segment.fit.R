
# Note: We should call this function whatever we end up naming the package
bulk.segment.fit = function(
  coverage.model.obj,
  eps = 10^-4,
  seed = NULL,
  truncated.models = FALSE,
  uniform.peak.threshold = 0.75,
  uniform.peak.stepsize = 5,
  remove.low.entropy = T,
  max.uniform = T,
  histogram.metric = c("jaccard", "intersection", "ks", "mse", "chisq")
){
  
  # Coverage.model.obj
  cov = coverage.model.obj$histogram.coverage
  histogram.ids = names(cov)
  
  results = list()
  # Running segmentation & distribution fitting
  # This part we can run in parallel
  for(i in histogram.ids){
    res = segment.fit.agnostic(
      x = cov[[i]],
      eps = eps,
      seed = seed,
      truncated.models = truncated.models,
      uniform.peak.threshold = uniform.peak.threshold,
      uniform.peak.stepsize = uniform.peak.stepsize,
      remove.low.entropy = remove.low.entropy,
      max.uniform = max.uniform,
      histogram.metric = histogram.metric
    )
    results[[i]] <- res
  }
  
  coverage.model.obj[['results']] <- results
  coverage.model.obj
}

# Check that this is compatible with Base 0 and Base 1 formatting
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

# Formatting results
format.results = function(
  coverage.model.obj
){
  
  # Coverage.model.obj
  gene.model = coverage.model.obj$gene.model
  histogram.ids = names(gene.model)
  histogram.bin.size = coverage.model.obj$histogram.bin.size
  results = coverage.model.obj$results
  
  formatted.results = GenomicRanges::GRanges()
  for(i in histogram.ids){
    x = gene.model[i]
    x = unlist(x)
    bins = GenomicRanges::tile(x = x, width = histogram.bin.size)
    bins = unlist(bins)
    res = results[i]
    
    # Extract ranges & collapse
    # This part is specific to the code we have now 
    # but can be easily changed to match the format of the results
    # Note: We can also add in the rest of the bed12 columns
    res = res[res$final == 1,]
    seg.starts = res$seg.start
    seg.ends = res$seg.end
    formatted.res = GenomicRanges::GRanges()
    for(j in 1:nrow(res)){
      f.res = GenomicRanges::reduce(bins[seg.starts[j]:seg.ends[j]])
      f.res = bed6tobed12(f.res)
      formatted.res = c(formatted.res, f.res)
    }
    S4Vectors::mcols(formatted.res)$name = paste0(i, "::", res$i)
    S4Vectors::mcols(formatted.res)[,c( "value", "metric", "dist", "params")] = res[,c( "value", "metric", "dist", "params")]
    formatted.results = c(formatted.results, formatted.res)
  }
  # Reordering the columns
  S4Vectors::mcols(formatted.results) = S4Vectors::mcols(formatted.results)[,c("name", "value", "metric", "dist", "params", "blockCount", "blockSizes", "blockStarts")]
  formatted.results
}