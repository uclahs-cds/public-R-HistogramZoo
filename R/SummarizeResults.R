results_columns = c(
  "region_id",
  "peak_id",
  "chr",
  "start",
  "end",
  "strand",
  "blockCount",
  "blockSizes",
  "blockStarts",
  "histogram_start",
  "histogram_end",
  "value",
  "metric",
  "dist",
  "params"
)

extract.stats = function(models, peak.id){
  mod = models$majority.vote
  df = data.frame(
    "peak_id" = peak.id,
    "histogram_start" = mod$seg.start,
    "histogram_end" = mod$seg.end,
    "value" = mod$value,
    "metric" = mod$metric,
    "dist" = mod$dist,
    "params" = dput.str(mod$par)
  )
  df
}

#' Formats results of SegmentAndFit
#'
#' @param result A Histogram object which have attributes 'models' and 'p' (i.e. the return object of running SegmentAndFit on a Histogram)
#'
#' @return TODO: Describe columns of dataframe
#' @export
SummarizeResults = function(result){
  UseMethod('SummarizeResults')
}

#' Title
#'
#' @param result
#'
#' @return
#' @export
#'
#' @examples
SummarizeResults.Histogram = function(
  result = NULL
){

  # Error checking
  stopifnot(inherits(result, "Histogram"))
  if(is.null(results$models)){
    stop("No models found. Please run SegmentAndFit first.")
  }

  # Attributes of the result
  models = result$models
  interval_start = result$interval_start
  interval_end = result$interval_end
  region_id = result$region_id

  # Generating a results table
  results_table = lapply(seq_along(models), function(i){
    extract.stats(models = models[[i]], peak.id = i)
  })
  results_table = do.call('rbind.data.frame', results_table)

  # Recovering coordinates for genes
  bins = IRanges::IRanges(start = interval_start, end = interval_end)
  bed.coords = apply(results_table, 1, function(peak){
    coords.gr = IRanges::reduce(bins[peak['histogram_start']:peak['histogram_end']])
    coords.gr = base1.to.base0(coords.gr)
    coords.gr = bed6tobed12(coords.gr)
    coords.gr
  })
  coords_table = do.call(c, bed.coords)
  coords_table = data.frame(coords_table)
  coords_table = subset(coords_table, select=-width)

  results_table = cbind.data.frame(results_table, coords_table)
  results_table['region_id'] <- region_id

  # Reorder columns
  results_table = results_table[,order(match(colnames(results_table), results_columns))]

  return(results_table)
}

#'
#' @param result TODO
#'
#'
#' @return TODO
#' @export
SummarizeResults.GenomicHistogram = function(
  result = NULL
){

  # Error checking
  stopifnot(inherits(result, "GenomicHistogram"))

  # Calling SummarizeResults.Histogram
  results_table = SummarizeResults.Histogram(result)

  # Add chromosome and strand for GenomicHistogram
  results_table['chr'] <- results$chr
  results_table['strand'] <- results$strand

  # Reorder columns
  results_table = results_table[,order(match(colnames(results_table), results_columns))]

  return(results_table)
}
