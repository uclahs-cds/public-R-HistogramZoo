results_columns <- c(
  "region_id",
  "segment_id",
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
  "sample_mean",
  "sample_var",
  "sample_sd",
  "sample_skew"
)

#' Extract stats from models
#'
#' @param model_list a list of models represented as a nested list, like the output of `fit_distributions`
#' @param model_name character, the name of the model from which the stats are to be extracted
#'
#' @return A list of the following data representing a single segment
#' \describe{
#'     \item{histogram_start}{The start index of the segment in the Histogram representation}
#'     \item{histogram_end}{The end index of the segment in the Histogram representation}
#'     \item{value}{The fitted value of the metric function}
#'     \item{metric}{The metric used to fit the distribution}
#'     \item{dist}{The distribution name}
#'     \item{params}{The parameters of the distribution}
#'}
extract_stats_from_models <- function(model_list, model_name = "consensus"){
  mod <- model_list[[model_name]]
  list(
    "histogram_start" = mod$seg_start,
    "histogram_end" = mod$seg_end,
    "value" = mod$value,
    "metric" = mod$metric,
    "dist" = mod$dist,
    "params" = mod$par
  )
}

#' Rescale parameter statistics based on original fit bin width
#'
#' @param stats a list of model statistics. Output from `extract_stats_from_models`
#' @param bin_width numeric, the size of the bin_width fit on the data
#'
#' @return A list of the following data representing a single segment
#' @export
scale_model_params <- function(stats, bin_width = 1){
  stopifnot(bin_width >= 1)
  if(stats$dist == 'norm') {
      stats$params$mean <- stats$params$mean * bin_width
      stats$params$sd <- stats$params$sd * bin_width
  } else if(stats$dist == 'gamma' || stats$dist == 'gamma_flip') {
    stats$params$rate <- stats$params$rate / bin_width
  }
  return(stats)
}

#' Represent a single segment as a set of intervals
#'
#' @param iranges an IRanges object with intervals representing intervals of a single segment
#'
#' @return a list with the following data representing a single segment
#' \describe{
#'     \item{start}{the interval start of the segment}
#'     \item{end}{the interval end of the segment}
#'     \item{interval_count}{The number of intervals in the segment - used for collapsing disjoint intervals}
#'     \item{interval_sizes}{The width of each interval}
#'     \item{interval_starts}{The start index of each interval}
#'}
extract_segments <- function(iranges){
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
#' @param result a Histogram object which have attributes 'models' and 'p' (i.e. the return object of running segment_and_fit on a Histogram)
#' @param model_name One of the metrics used to fit models (e.g. Jaccard) and "consensus" if more than one metric was used to specify which model params to extract
#'
#' @return a data.frame with all or a subset of the following columns summarizing the results of the fit
#' \describe{
#'     \item{region_id}{character string denoting the region_id of the Histogram}
#'     \item{segment_id}{an integer id identifying the ordinal segment of the Histogram segmentation}
#'     \item{chr}{an optional column denoting the chromosome of a GenomicHistogram object}
#'     \item{start}{the interval start of the segment}
#'     \item{end}{the interval end of the segment}
#'     \item{strand}{an optional column denoting the strand of a GenomicHistogram object}
#'     \item{interval_count}{the number of intervals in the segment - used for collapsing disjoint intervals}
#'     \item{interval_sizes}{the width of each interval}
#'     \item{interval_starts}{the start index of each interval}
#'     \item{histogram_start}{the start index of the segment in the Histogram representation}
#'     \item{histogram_end}{the end index of the segment in the Histogram representation}
#'     \item{value}{the fitted value of the metric function}
#'     \item{metric}{the metric used to fit the distribution}
#'     \item{dist}{the distribution name}
#'     \item{sample_mean}{weighted mean of the histogram}
#'     \item{sample_var}{weighted variance of the histogram}
#'     \item{sample_sd}{weighted standard deviation of the histogram}
#'     \item{sample_skew}{weighted skew of the histogram}
#'     \item{dist_param[0-9]}{the values of distribution parameters}
#'     \item{dist_param_name[0-9]}{the matching names of distribution parameters}
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
  model_name <- match.arg(model_name, c("consensus", result$metric))

  # Attributes of the result
  models <- result$models
  interval_start <- result$interval_start
  interval_end <- result$interval_end
  region_id <- result$region_id

  # Summarizing results
  max_param_length <- max(unlist(lapply(models, function(m) {
    length(m[[model_name]]$par)
  })))
  if(max_param_length > 0) {
    param_names <- paste0('dist_param', seq_len(max_param_length))
    params <- rep(NA, max_param_length)
    names(params) <- param_names
  }

  bins <- IRanges::IRanges(start = interval_start, end = interval_end)
  results_table <- lapply(seq_along(models), function(i){

    # Model parameters
    stats <- extract_stats_from_models(model_list = models[[i]], model_name = model_name)
    # if(result$bin_width > 1) {
    #   stats <- scale_model_params(stats, bin_width = result$bin_width)
    # }

    # Bin coordinates
    coords <- IRanges::reduce(bins[stats[['histogram_start']]:stats[['histogram_end']]])
    coords <- extract_segments(coords)

    # Empirical moment estimation from histogram
    moment_estimation <- list(
      "sample_mean" = weighted.mean(result[stats[['histogram_start']]:stats[['histogram_end']]]),
      "sample_var" = weighted.var(result[stats[['histogram_start']]:stats[['histogram_end']]]),
      "sample_sd" = weighted.sd(result[stats[['histogram_start']]:stats[['histogram_end']]]),
      "sample_skew" = weighted.skewness(result[stats[['histogram_start']]:stats[['histogram_end']]])
    )

    # Parameters
    if (max_param_length > 0) {
      # Convert parameters into common naming scheme rather than
      # actual parameter names like mean, sd, rate, shape, etc
      new_params <- params
      new_param_names <- rep(NA, max_param_length)
      n_params <- length(stats$params)
      if(n_params > 1) {
        old_params <- stats$params
        new_params[1:n_params] <- old_params
        new_param_names[1:n_params] <- names(stats$params)
      }
      new_params <- c(new_params, new_param_names)
      names(new_params) <- c(param_names, paste0(param_names, '_name'))

      stats$params <- NULL

      data.frame(stats[lengths(stats) > 0], as.list(new_params), coords, moment_estimation, segment_id = i)
    } else {
      data.frame(stats[lengths(stats) > 0], coords, moment_estimation, segment_id = i)
    }
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
