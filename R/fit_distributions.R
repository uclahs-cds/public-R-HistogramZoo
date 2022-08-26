#' Fit a uniform distribution to a histogram
#'
#' @param x numeric vector representing the density of a histogram
#' @param metric a subset of `jaccard`, `intersection`, `ks`, `mse`, `chisq`
#' 
#' @return A list with the following data
#' \describe{
#'     \item{par}{A character string denoting the region_id of the Histogram}
#'     \item{dist}{The distribution name}
#'     \item{metric}{The metric used to fit the distribution}
#'     \item{value}{The fitted value of the metric function}
#'     \item{dens}{A function that returns the density of the fitted distribution}
#' }
fit_uniform <- function(x, metric){
  
  N <- sum(x)
  bin <- 1:length(x)
  p_unif <- generate_uniform_distribution(x)
  h_unif <- x / sum(x)
  metric_func <- get(paste('histogram', metric, sep = "."))
  
  m <- metric_func(h_unif, p_unif)

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
#' @return A nested list where each sublist represents a model with the following data
#' \describe{
#'     \item{par}{A character string denoting the region_id of the Histogram}
#'     \item{dist}{The distribution name}
#'     \item{metric}{The metric used to fit the distribution}
#'     \item{value}{The fitted value of the metric function}
#'     \item{dens}{A function that returns the density of the fitted distribution}
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
  .hist.optim <- function(params, .dist = c("norm", "gamma")) {
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
    res <- metric_func(histogram_data, dens)
    if(is.na(res) || res == -Inf || res == Inf) {
      res <- Inf
    }
    res
  }


  # Initializing Todo & Results
  todo <- expand.grid(metric, distributions, stringsAsFactors = F)
  rtn <- list()
  for(i in 1:nrow(todo)){

    # Extracting Metric & Distribution
    met <- todo[i,1]
    dist <- todo[i,2]
    tag <- paste0(met, ".", dist)

    # Get one of the metrics from histogram.distances
    metric_func <- get(paste('histogram', met, sep = "."))

    # Uniform Distribution
    if(dist == "unif") {
      # Add uniform distribution
      unif_dens <- 1 / (max(bin) - min(bin))
      rtn[[tag]] <- list(
        "par" = NULL,
        "dist" = "unif",
        "metric" = met,
        "value" = correct_fitted_value(met, metric_func(histogram_data, rep(unif_dens * N, L))),
        "dens" = function(x = NULL, mpar = NULL, scale = TRUE) {
          if(missing(x)) {
            x <- bin
          }
          res <- ifelse(x >= min(bin) & x <= max(bin), unif_dens, 0)
          if(scale) res * N
          else res
        }
      )
    }
    # Normal Distribution
    if(dist == "norm") {
      norm_optim <- DEoptim::DEoptim(fn = .hist.optim,
                                    .dist = "norm",
                                    lower = c(min(bin), 0.001),
                                    upper = c(max(bin), (max(bin) - min(bin)) * 0.5),
                                    control = list(
                                    trace = FALSE, # Do not print results
                                    itermax = 500, # Iterations
                                    VTR = 10^-2 # At 1 %, stop optimizing
                                    ))
      names(norm_optim$optim$bestmem) <- c("mean", "sd")
      norm_par <- as.list(norm_optim$optim$bestmem)
      if(truncated) {
        norm_par$a <-  min(bin) - 1e-10
        norm_par$b <- max(bin) + 1e-10
      }
      rtn[[tag]] <- list(
        "par" = norm_par,
        "dist" = "norm",
        "metric" = met,
        "value" = correct_fitted_value(met, norm_optim$optim$bestval),
        "dens" = function(x = NULL, mpar = NULL, scale = TRUE) {
          if(missing(x)) {
            x <- bin
          }
          args <- c(list(x = x), as.list(mpar))
          res <- do.call("dtnorm", args)
          if(scale) res * N
          else res
        }
      )
    }

    # Gamma Distribution
    if(dist == "gamma") {
      gamma_optim <- DEoptim::DEoptim(fn = .hist.optim,
                                     .dist = "gamma",
                                     lower = c(0.001, 0.001),
                                     upper = c(L, L),
                                     control = list(
                                     trace = FALSE, # Do not print results
                                     itermax = 500, # Iterations
                                     VTR = 10^-2 # At 1 %, stop optimizing
                                     ))
      names(gamma_optim$optim$bestmem) <- c("shape", "rate")
      gamma_par <- as.list(gamma_optim$optim$bestmem)
      if(truncated) {
        gamma_par$a <- min(bin) - 1e-10
        gamma_par$b <- max(bin) + 1e-10
      }
      rtn[[tag]] <- list(
        "par" = gamma_par,
        "dist" = "gamma",
        "metric" = met,
        "value" = correct_fitted_value(met, gamma_optim$optim$bestval),
        "dens" = function(x = NULL, mpar = NULL, scale = TRUE) {
          if(missing(x)) {
            x <- bin
          }
          args <- c(list(x = x), as.list(mpar))
          res <- do.call("dtgamma", args)
          if(scale) res * N
          else res
        }
      )
  }
  }

  return(rtn)
}
