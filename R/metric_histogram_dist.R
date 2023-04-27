#' Computes a given histogram distance metric to a distribution
#'
#' Generally this function is used to estimate the parameters of a distribution
#' by minimizing the given metric. The first argument is a vector of parameters
#' that should match up with the arguments for `dist`. For example, \code{params = c(0, 1)} with
#' \code{dist = 'norm'} refers to mean 0 and sd 1. This is so that we can use
#' \code{metric.histogram.dist} in [`DEoptim::DEoptim`] which requires the parameters to
#' optimize to be the first argument.
#'
#' @param params vector of parameters that should match up with the arguments for `dist`.
#' @param x counts of the histogram
#' @param interval_start left endpoint of the histogram
#' @param interval_end right endpoint of the histogram
#' @param dist one of `norm`, `gamma` or `gamma_flip`
#' @param metric_func one of
#' `r paste0(ls(pattern = '^histogram[.].*', envir = getNamespace('HistogramZoo')), collapse = ', ')`
#' @param truncated Should distributions be truncated to the end points of the histogram?
#'
#' @return scalar results from comparing the distribution, with given parameters to observed histogram
#' @export
#'
#' @examples
#' x.hist <- observations_to_histogram(rnorm(1000))
#'
#' args <- list(
#'   x = x.hist$histogram_data,
#'   interval_start = x.hist$interval_start,
#'   interval_end = x.hist$interval_end,
#'   dist = 'norm',
#'   metric_func = 'histogram.intersection'
#'   )
#'
#' do.call(
#'   metric.histogram.dist,
#'   c(args, list(params = c(0, 1)))
#'   )
#'
#' do.call(
#'   metric.histogram.dist,
#'   c(args, list(params = c(1, 1.5)))
#'   )
metric.histogram.dist <- function(
    params,
    x,
    interval_start,
    interval_end,
    dist = c("norm", "gamma", "gamma_flip"),
    metric_func = ls(pattern = '^histogram[.].*', envir = getNamespace('HistogramZoo')),
    truncated = FALSE
  ) {
  possible.metrics <- ls(pattern = '^histogram[.].*', envir = getNamespace('HistogramZoo'))
  if (missing(metric_func)) {
    stop(
      'Provide a metric function for argument `metric_func`: ', paste0(possible.metrics, collapse = ', ')
      )
  }

  area <- sum(x*(interval_end - interval_start))
  interval_midpoint <- (interval_end + interval_start) / 2

  if (is.character(metric_func)) metric_func <- get(metric_func)

  # Compute the expected counts for the given parameters
  args <- c(list(x = interval_midpoint), params)
  if(truncated && dist == "normal") {
    args$a <- head(interval_start, 1) - 1e-10
    args$b <- tail(interval_end, 1) + 1e-10
  }
  if(truncated && dist %in% c("gamma_flip", "gamma")){
    args$a <- 0
    args$b <- L
  }
  if(dist == "gamma_flip"){
    args$offset <- tail(interval_end, 1) + 1e-10
  }
  if(dist == "gamma"){
    args$shift <- head(interval_start, 1) - 1e-10
  }
  trunc.letter <- if(truncated) "t" else ""
  dens <- tryCatch({
    do.call(paste0("d", trunc.letter, dist), args) * area
  }, error = function(err) {
    rep(0, length(x))
  })
  dens[is.na(dens)] <- 0
  res <- metric_func(x, dens)
  if(is.na(res) || res == -Inf || res == Inf) {
    res <- Inf
  }
  res
}
