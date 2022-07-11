
extract.stats = function(models, gene = NULL, peak.id){
  mod = models$majority.vote
  df = data.frame(
    "peak.id" = peak.id,
    "start" = mod$seg.start,
    "end" = mod$seg.end,
    "value" = mod$value,
    "metric" = mod$metric,
    "dist" = mod$dist,
    "params" = dput.str(mod$par)
  )
  if(!is.null(gene)) {df = cbind("gene" = gene, df)}
  df
}

# Method Dispatch
SummarizeResults = function(result){
  UseMethod('SummarizeResults', result)
}

#' Formats results of segment.fit.agnostic
#'
#' @param result TODO
#'
#' @return TODO
#' @export
SummarizeResults.Histogram = function(
  result = NULL
){

  # Error checking
  stopifnot(inherits(result, "Histogram"))
  if(is.null(attr(result, "models")) | is.null(attr(result, "p"))){
    stop("No computed results. Please run SegmentAndFit first.")
  }

  res = attr(result, "models")
  results.tbl = lapply(seq_along(res), function(i){
    extract.stats(
      models = res[[i]],
      gene = NULL,
      peak.id = i)
  })
  results.tbl = do.call('rbind.data.frame', results.tbl)
  results.tbl
}

#' Formats results of bulk.segment.fit
#'
#' @param result TODO
#' @param output.format TODO
#'
#'
#' @return TODO
#' @export
SummarizeResults.GenomicHistogram = function(
  result = NULL
){

  # Error checking
  stopifnot(inherits(result, "GenomicHistogram"))
  if(is.null(attr(result, "models")) | is.null(attr(result, "p"))){
    stop("No computed results. Please run SegmentAndFit first.")
  }

  # Recovering output stats
  res = attr(result, "models")
  results.tbl = lapply(seq_along(res), function(i){
    extract.stats(
      models = res[[i]],
      gene = NULL,
      peak.id = i)
  })
  results.tbl = do.call('rbind.data.frame', results.tbl)

  # Recovering coordinates for genes
  interval_start = attr(result, "interval_start")
  interval_end = attr(result, "interval_end")
  bins = IRanges::IRanges(start = interval_start, end = interval_end)
  bed.coords = apply(results.tbl, 1, function(peak){
    coords.gr = IRanges::reduce(bins[peak['start']:peak['end']])
    coords.gr = base1.to.base0(coords.gr)
    coords.gr = bed6tobed12(coords.gr)
    coords.gr
  })
  coords.tbl = do.call(c, bed.coords)
  coords.tbl = data.frame(coords.tbl)
  coords.tbl = subset(coords.tbl, select=-width)

  # Technically this is the only part that might be unique to GenomicHistogram
  coords.tbl['chr'] <- attr(results, "chr")
  coords.tbl['strand'] <- attr(results, "strand")

  # TODO: Potentially rewrite the bottom so that this is more efficient
  results.tbl = subset(results.tbl, select=-c(start, end))
  results.tbl = cbind.data.frame(results.tbl, coords.tbl)
  
  results.tbl
}
