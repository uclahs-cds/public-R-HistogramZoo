#' Fit a uniform distribution to a histogram
#'
#' @param x numeric vector representing the density of a histogram
#' @param metric a subset of `jaccard`, `intersection`, `ks`, `mse`, `chisq`
#'
#' @return a list with the following data
#' \describe{
#'     \item{par}{a character string denoting the region_id of the Histogram}
#'     \item{dist}{the distribution name}
#'     \item{metric}{the metric used to fit the distribution}
#'     \item{value}{the fitted value of the metric function}
#'     \item{dens}{a function that returns the density of the fitted distribution}
#' }
fit_uniform <- function(x, metric){

  N <- sum(x)
  bin <- 1:length(x)
  p_unif <- generate_uniform_distribution(x)
  # h_unif <- x / sum(x)
  metric_func <- get(paste('histogram', metric, sep = "."))

  # m <- metric_func(h_unif, p_unif)
  m <- metric_func(x, p_unif*N)

  return(
    list(
      "par" = NULL,
      "dist" = "unif",
      "metric" = metric,
      "value" = correct_fitted_value(metric, m),
      "dens" = function(x = NULL, mpar = NULL, scale = TRUE) {
        if(missing(x)) {
          x <- bin
        }
        res <- ifelse(x >= min(bin) & x <= max(bin), p_unif[1], 0)
        if(scale) res * N
        else res
      }
    )
  )

}

#' Fit the model parameters by optimizing a histogram metric
#'
#' @param histogram_data numeric vector, representing data to be fit
#' @param metric a subset of `jaccard`, `intersection`, `ks`, `mse`, `chisq`
#' indicating metrics to use for fit optimization
#' @param truncated logical, whether to fit truncated distributions
#' @param distributions character vector indicating distributions,
#' subset of `norm`, `gamma`, `unif`
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
#' @importFrom DEoptim DEoptim
fit_distributions <- function(
    histogram_data,
    metric = c("jaccard", "intersection", "ks", "mse", "chisq"),
    truncated = FALSE,
    distributions = c("norm", "gamma", "unif")) {

  # Matching arguments
  # TODO: consider checking minimum length of histogram_data or
  # check if Histogram object
  if(!is.numeric(histogram_data)){
    stop('histogram_data has to be a numeric vector')
  }
  if(!is.logical(truncated) | length(truncated) != 1){
    stop("truncated has to be a logical of length 1")
  }
  metric <- match.arg(metric, several.ok = TRUE)
  distributions <- match.arg(distributions, several.ok = TRUE)

  # Initializing Data
  L <- length(histogram_data)
  N <- sum(histogram_data)
  bin <- 1:L

  # Optimization Function
  .hist.optim <- function(params, .dist = c("norm", "gamma"), .metric_func) {
    # Compute the expected counts for the given parameters
    args <- c(list(x = bin), params)
    if(truncated) {
      args$a <-  min(bin) - 1e-10
      args$b <- max(bin) + 1e-10
    }
    trunc.letter <- if(truncated) "t" else ""
    dens <- tryCatch({
      do.call(paste0("d", trunc.letter, .dist), args) * N
    }, error = function(err) {
      rep(0, length(bin))
    })
    dens[is.na(dens)] <- 0
    res <- .metric_func(histogram_data, dens)
    if(is.na(res) || res == -Inf || res == Inf) {
      res <- Inf
    }
    res
  }

  # Fitting uniform distributions
  unif_fit <- list()
  if("unif" %in% distributions){
    distributions <- setdiff(distributions, "unif")
    unif_fit <- lapply(metric, function(met) fit_uniform(histogram_data, met))
  }

  rtn <- lapply(distributions, function(distr) {
    lapply(metric, function(met){

      # Get one of the metrics from histogram.distances
      metric_func <- get(paste('histogram', met, sep = "."))

      # Setting boundaries & parameter names
      if(distr == "norm"){
        lower <- c(min(bin), 0.001)
        upper <- c(max(bin), (max(bin) - min(bin)) * 0.5)
        names_par <- c("mean", "sd")
      } else if (distr == "gamma"){
        lower <- c(0.001, 0.001)
        upper <- c(L, L)
        names_par <- c("shape", "rate")
      }

      # Fitting Data
      dist_optim <- DEoptim::DEoptim(
        fn = .hist.optim,
        .dist = distr,
        .metric_func = metric_func,
        lower = lower,
        upper = upper,
        control = list(
          trace = FALSE, # Do not print results
          itermax = 500, # Iterations
          VTR = 10^-2 # At 1 %, stop optimizing
        )
      )

      # Extracting parameters
      names(dist_optim$optim$bestmem) <- names_par
      dist_par <- as.list(dist_optim$optim$bestmem)

      # Adjusting for truncated distributions
      if(truncated) {
        dist_par$a <-  min(bin) - 1e-10
        dist_par$b <- max(bin) + 1e-10
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
              x <- bin
            }
            args <- c(list(x = x), as.list(mpar))
            res <- do.call(paste0("dt", distr), args)
            if(scale) res * N
            else res
          }
        )
      )
      ##
    })
  })

  return(
    c( unlist(rtn, recursive = F), unif_fit )
  )
}
