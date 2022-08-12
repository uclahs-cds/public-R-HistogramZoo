results_columns <- c(
  "region_id",
  "peak_id",
  "chr",
  "start",
  "end",
  "strand",
  "interval_count",
  "interval_sizes",
  "interval_starts",
  "histogram_start",
  "histogram_end",
  "value",
  "metric",
  "dist",
  "params"
)

extract_stats_from_models <- function(model_list, model_name = "consensus"){
  mod <- model_list[[model_name]]
  list(
    "histogram_start" = mod$seg_start,
    "histogram_end" = mod$seg_end,
    "value" = mod$value,
    "metric" = mod$metric,
    "dist" = mod$dist,
    "params" = dput_str(mod$par)
  )
}

extract_peak_segments <- function(iranges){
  iranges_range <- range(iranges)
  list(
    "start" = IRanges::start(iranges_range),
    "end" = IRanges::end(iranges_range),
    "interval_count" = length(iranges),
    "interval_sizes" = paste0(IRanges::width(iranges), ",", collapse = ""),
    "interval_starts" = paste0(IRanges::start(iranges) - IRanges::start(iranges_range) + 1, ",", collapse = "")
  )
}

#' Formats results of segment_and_fit
#'
#' @param result A Histogram object which have attributes 'models' and 'p' (i.e. the return object of running segment_and_fit on a Histogram)
#' @param model_name One of the metrics used to fit models (e.g. Jaccard) and "consensus" if more than one metric was used to specify which model params to extract
#'
#' @return A data.frame with the following columns summarizing the results of the fit
#' \describe{
#'     \item{region_id}{character string denoting the region_id of the Histogram}
#'     \item{peak_id}{an integer id identifiying the ordinal segment of the Histogram segmentation}
#'     \item{chr}{an optional column denoting the chromosome of a GenomicHistogram object}
#'     \item{start}{the interval start of the segment}
#'     \item{end}{the interval end of the segment}
#'     \item{strand}{an optional column denoting the strand of a GenomicHistogram object}
#'     \item{interval_count}{}
#'     \item{interval_sizes}{}
#'     \item{interval_starts}{}
#'     \item{histogram_start}{}
#'     \item{histogram_end}{}
#'     \item{value}{}
#'     \item{metric}{}
#'     \item{dist}{}
#'     \item{params}{}
#' }
#'
#' @export
summarize_results <- function(result, model_name){
  UseMethod('summarize_results')
}


#' @export
summarize_results.Histogram <- function(
  result = NULL,
  model_name = "consensus"
){

  # Error checking
  stopifnot(inherits(result, "HistogramFit"))
  stopifnot(inherits(result, "Histogram"))
  model_name <- match.arg(model_name, c("consensus", result$histogram_metric))

  # Attributes of the result
  models <- result$models
  interval_start <- result$interval_start
  interval_end <- result$interval_end
  region_id <- result$region_id

  # Summarizing results
  bins <- IRanges::IRanges(start = interval_start, end = interval_end)
  results_table <- lapply(seq_along(models), function(i){
    stats <- extract_stats_from_models(model_list = models[[i]], model_name = model_name)
    coords <- IRanges::reduce(bins[stats[['histogram_start']]:stats[['histogram_end']]])
    coords <- extract_peak_segments(coords)
    c(stats, coords, list('peak_id' = i))
  })
  results_table <- do.call('rbind.data.frame', results_table)
  results_table['region_id'] <- region_id

  # Reorder columns
  results_table <- results_table[,order(match(colnames(results_table), results_columns))]

  return(results_table)
}

#' @export
summarize_results.GenomicHistogram = function(
  result = NULL,
  model_name = "consensus"
){

  # Error checking
  stopifnot(inherits(result, "HistogramFit"))
  stopifnot(inherits(result, "GenomicHistogram"))

  # Calling summarize_results.Histogram
  results_table <- summarize_results.Histogram(result, model_name)

  # Add chromosome and strand for GenomicHistogram
  results_table['chr'] <- result$chr
  results_table['strand'] <- result$strand

  # Reorder columns
  results_table <- results_table[,order(match(colnames(results_table), results_columns))]

  return(results_table)
}
