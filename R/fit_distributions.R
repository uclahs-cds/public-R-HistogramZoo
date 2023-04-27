
#' Fit the model parameters by optimizing a histogram metric
#'
#' @param x numeric vector, representing data to be fit
#' @param metric a subset of `mle`, `jaccard`, `intersection`, `ks`, `mse`, `chisq`
#' indicating metrics to use for fit optimization
#' @param truncated logical, whether to fit truncated distributions
#' @param distributions character vector indicating distributions,
#' subset of `norm`, `gamma`, `gamma_flip` and `unif`.
#'
#' @export
#'
#' @return a nested list where each sublist represents a model with the following data
#' \describe{
#'     \item{par}{a character string denoting the region_id of the Histogram}
#'     \item{dist}{the distribution name}
#'     \item{metric}{the metric used to fit the distribution}
#'     \item{value}{the fitted value of the metric function}
#'     \item{dens}{a function that returns the density of the fitted distribution}
#' }
#'
#' @rdname fit_distributions
#' @importFrom DEoptim DEoptim
fit_distributions <- function(
    x,
    metric = c("mle", "jaccard", "intersection", "ks", "mse", "chisq"),
    truncated = FALSE,
    distributions = c("norm", "gamma", "gamma_flip", "unif")
) {
  UseMethod("fit_distributions")
}

#' @rdname fit_distributions
#' @exportS3Method fit_distributions numeric
fit_distributions.numeric <- function(
    x,
    metric = c("mle", "jaccard", "intersection", "ks", "mse", "chisq"),
    truncated = FALSE,
    distributions = c("norm", "gamma", "gamma_flip", "unif")
) {

  # Assume a bin width of 1 for a numeric vector
  L <- length(x)
  interval_midpoint <- seq(from = 1, to = L, by = 1)
  fit_distributions_helper(
    x = x,
    interval_start = interval_midpoint - 0.5,
    interval_end = interval_midpoint + 0.5,
    interval_midpoint = interval_midpoint,
    metric = metric,
    truncated = truncated,
    distributions = distributions
  )

}

#' @rdname fit_distributions
#' @exportS3Method fit_distributions table
fit_distributions.table <- fit_distributions.numeric

#' @rdname fit_distributions
#' @exportS3Method fit_distributions GenomicHistogram
fit_distributions.GenomicHistogram <- function(
    x,
    metric = c("mle", "jaccard", "intersection", "ks", "mse", "chisq"),
    truncated = FALSE,
    distributions = c("mle", "norm", "gamma", "gamma_flip", "unif")
) {

  # Base 1 to a base 0 system for GenomicHistogram
  fit_distributions_helper(
    x = x$histogram_data,
    interval_start = x$consecutive_start - 0.5,
    interval_end = x$consecutive_end + 0.5,
    interval_midpoint = find_midpoint(x),
    metric = metric,
    truncated = truncated,
    distributions = distributions
  )

}

#' @rdname fit_distributions
#' @exportS3Method fit_distributions Histogram
fit_distributions.Histogram <- function(
    x,
    metric = c("mle", "jaccard", "intersection", "ks", "mse", "chisq"),
    truncated = FALSE,
    distributions = c("norm", "gamma", "gamma_flip", "unif")
) {

  # Interval start and end from Histogram
  fit_distributions_helper(
    x = x$histogram_data,
    interval_start = x$interval_start,
    interval_end = x$interval_end,
    interval_midpoint = find_midpoint(x),
    metric = metric,
    truncated = truncated,
    distributions = distributions
  )

}

#' A helper function for fit_distributions
#'
#' @param x numeric vector, representing data to be fit
#' @param interval_start starting coordinates for the bins
#' @param interval_end ending coordinates for the bins
#' @param interval_midpoint bin midpoints
#' @param metric a subset of `mle`, `jaccard`, `intersection`, `ks`, `mse`, `chisq`
#' indicating metrics to use for fit optimization
#' @param truncated logical, whether to fit truncated distributions
#' @param distributions character vector indicating distributions,
#' subset of `norm`, `gamma`, `gamma_flip` and `unif`.
#'
#' @return a nested list where each sublist represents a model with the following data
#' \describe{
#'     \item{par}{a character string denoting the region_id of the Histogram}
#'     \item{dist}{the distribution name}
#'     \item{metric}{the metric used to fit the distribution}
#'     \item{value}{the fitted value of the metric function}
#'     \item{dens}{a function that returns the density of the fitted distribution}
#' }
#'
#' @importFrom DEoptim DEoptim
fit_distributions_helper <- function(
    x,
    interval_start,
    interval_end,
    interval_midpoint,
    metric = c("mle", "jaccard", "intersection", "ks", "mse", "chisq"),
    truncated = FALSE,
    distributions = c("norm", "gamma", "gamma_flip", "unif")
) {

  # Matching arguments
  # NOTE: not checking validity of interval start/end/midpoint
  # because computed internally
  if(!is.logical(truncated) | length(truncated) != 1){
    stop("truncated has to be a logical of length 1")
  }
  metric <- match.arg(
    metric,
    several.ok = TRUE,
    choices = c("mle", "jaccard", "intersection", "ks", "mse", "chisq")
  )
  distributions <- match.arg(distributions, several.ok = TRUE)

  # Setting vars for parameter estimation
  L <- tail(interval_end, 1) - head(interval_start, 1)
  # Area = sum(density * bin_width)
  area <- sum(x*(interval_end - interval_start))

  # Fitting uniform distributions
  unif_fit <- list()
  if("unif" %in% distributions){
    distributions <- setdiff(distributions, "unif")
    unif_fit <- lapply(metric, function(met) fit_uniform_helper(
      x,
      interval_start,
      interval_end,
      interval_midpoint,
      met
    ))
  }

  rtn <- lapply(distributions, function(distr) {
    lapply(metric, function(met){
      # Setting boundaries & parameter names
      if(distr == "norm"){
        lower <- c(head(interval_start, 1), 0.001)
        upper <- c(tail(interval_end, 1), (tail(interval_end, 1) - head(interval_start, 1)) * 0.5)
        names_par <- c("mean", "sd")
      } else if (distr %in% c("gamma", "gamma_flip")){
        lower <- c(0.001, 0.001)
        upper <- c(L, L)
        names_par <- c("shape", "rate")
      }

      control_args <- list(
        trace = FALSE,
        itermax = 500,
        steptol = 50
        )

      if (met == 'mle') {
        tdistr <- distr
        mle.options.args <- list(
          fn = bin_log_likelihood,
          counts = x,
          bin_lower = interval_start,
          bin_upper = interval_end,
          lower = lower,
          upper = upper,
          control = control_args
          )

        if (distr == "gamma") {
          mle.options.args$shift <- head(interval_start, 1) - 1e-10
        } else if (distr == "gamma_flip") {
          mle.options.args$offset <- tail(interval_end, 1) + 1e-10
        }

        if (truncated) {
          tdistr <- paste0('t', tdistr)
          if(distr == "normal"){
            mle.options.args$a <- head(interval_start, 1) - 1e-10
            mle.options.args$b <- tail(interval_end, 1) + 1e-10
          } else if (distr %in% c("gamma", "gamma_flip")){
            mle.options.args$a <- 0
            mle.options.args$b <- L
          }

        }

        mle.options.args$cdf = get(paste0('p', tdistr))

        dist_optim <- do.call(DEoptim::DEoptim, mle.options.args)
      } else {
        # Get one of the metrics from histogram.distances
        metric_func <- get(paste('histogram', met, sep = "."))

        # Fitting Data
        dist_optim <- DEoptim::DEoptim(
          fn = metric.histogram.dist,
          x = x,
          interval_start = interval_start,
          interval_end = interval_end,
          dist = distr,
          metric_func = metric_func,
          lower = lower,
          upper = upper,
          control = control_args
        )
      }

      # Extracting parameters
      names(dist_optim$optim$bestmem) <- names_par
      dist_par <- as.list(dist_optim$optim$bestmem)

      # Adjusting for truncated distributions
      trunc.letter = if(truncated) "t" else ""
      if(truncated && distr == "normal") {
        dist_par$a <- head(interval_start, 1)
        dist_par$b <- tail(interval_end, 1)
      }
      # Adjusting for shift & offset
      if(distr == "gamma_flip"){
        dist_par$offset <- tail(interval_end, 1) + 1e-10
      }
      if(distr == "gamma"){
        dist_par$shift <- head(interval_start, 1) - 1e-10
      }
      if(truncated && distr %in% c("gamma", "gamma_flip")){
        dist_par$a <- 0
        dist_par$b <- L
      }

      # Return model
      return(
        list(
          "par" = dist_par,
          "dist" = distr,
          "metric" = met,
          "value" = correct_fitted_value(met, dist_optim$optim$bestval),
          "dens" = function(x = NULL, mpar = NULL, scale = TRUE) {
            if(missing(x)) {
              x <- interval_midpoint
            }
            args <- c(list(x = x), as.list(mpar))
            res <- do.call(paste0("d", trunc.letter, distr), args)
            if(scale) res * area
            else res
          }
        )
      )
    })
  })

  return(
    c(unlist(rtn, recursive = F), unif_fit)
  )
}
