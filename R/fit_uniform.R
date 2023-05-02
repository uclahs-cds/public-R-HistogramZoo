#' Fit a uniform distribution to a histogram
#'
#' @param x numeric vector representing the density of a histogram
#' @param metric one of `mle`, `jaccard`, `intersection`, `ks`, `mse`, `chisq`
#'
#' @return a list with the following data
#' \describe{
#'     \item{par}{a character string denoting the region_id of the Histogram}
#'     \item{dist}{the distribution name}
#'     \item{metric}{the metric used to fit the distribution}
#'     \item{value}{the fitted value of the metric function}
#'     \item{dens}{a function that returns the density of the fitted distribution}
#' }
#'
#' @rdname fit_uniform
#' @export
fit_uniform <- function(x, metric=c('mle', 'jaccard', 'intersection', 'ks', 'mse', 'chisq')) {
  UseMethod('fit_uniform')
}

#' @rdname fit_uniform
#' @exportS3Method fit_uniform numeric
fit_uniform.numeric <- function(x, metric=c('mle', 'jaccard', 'intersection', 'ks', 'mse', 'chisq')) {

  interval_midpoint <- seq(from = 1, to = length(x), by = 1)
  fit_uniform_helper(
    x = x,
    interval_start = interval_midpoint - 0.5,
    interval_end = interval_midpoint + 0.5,
    interval_midpoint = interval_midpoint,
    metric = metric
  )

}

#' @rdname fit_uniform
#' @exportS3Method  fit_uniform table
fit_uniform.table <- fit_uniform.numeric

#' @rdname fit_uniform
#' @exportS3Method fit_uniform Histogram
fit_uniform.Histogram <- function(x, metric=c('mle', 'jaccard', 'intersection', 'ks', 'mse', 'chisq')) {

  fit_uniform_helper(
    x = x$histogram_data,
    interval_start = x$interval_start,
    interval_end = x$interval_end,
    interval_midpoint = find_midpoint(x),
    metric = metric
  )

}

#' @rdname fit_uniform
#' @exportS3Method fit_uniform GenomicHistogram
fit_uniform.GenomicHistogram <- function(x, metric=c('mle', 'jaccard', 'intersection', 'ks', 'mse', 'chisq')) {

  fit_uniform_helper(
    x = x$histogram_data,
    interval_start = x$consecutive_start - 0.5,
    interval_end = x$consecutive_end + 0.5,
    interval_midpoint = find_midpoint(x),
    metric = metric
  )

}

#' A helper function for the S3 method fit_uniform
#'
#' @param x numeric vector representing the density of a histogram
#' @param interval_start starting coordinates for the bins
#' @param interval_end ending coordinates for the bins
#' @param interval_midpoint bin midpoints
#' @param metric one of `mle`, `jaccard`, `intersection`, `ks`, `mse`, `chisq`
#'
#' @return a list with the following data
#' \describe{
#'     \item{par}{a character string denoting the region_id of the Histogram}
#'     \item{dist}{the distribution name}
#'     \item{metric}{the metric used to fit the distribution}
#'     \item{value}{the fitted value of the metric function}
#'     \item{dens}{a function that returns the density of the fitted distribution}
#' }
fit_uniform_helper <- function(
    x,
    interval_start,
    interval_end,
    interval_midpoint,
    metric=c('mle', 'jaccard', 'intersection', 'ks', 'mse', 'chisq')
){

  metric <- match.arg(metric)

  # Initializing
  L <- tail(interval_end, 1) - head(interval_start, 1)
  area <- sum(x*(interval_end - interval_start))

  if(metric == 'mle') {
    # Negative log-likelihood
    value <- (- uniform_mle(
      x = x,
      x.start = interval_start,
      x.end = interval_end,
      a = head(interval_start, 1),
      b = tail(interval_end, 1),
      inclusive = TRUE,
      log = TRUE
    ))
  } else {
    p_unif <- rep(1/L, length(x))
    metric_func <- get(paste('histogram', metric, sep = "."))
    m <- metric_func(x, p_unif*area)
    value <- correct_fitted_value(metric, m)
  }

  return(
    new_ModelFit(
      par = NULL,
      dist = "unif",
      metric = metric,
      value = value,
      dens = function(x = NULL, mpar = NULL, scale = TRUE) {
        if(missing(x)) {
          x <- interval_midpoint
        }
        res <- ifelse(x >= head(interval_start, 1) & x <= tail(interval_end, 1), 1/L, 0)
        if(scale) res * area
        else res
      }
    )
  )

}
