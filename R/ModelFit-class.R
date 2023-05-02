
# constructor
#' Constructs a new HZModelFit object
#'
#' @param par list of named parameters
#' @param dist character name of distribution
#' @param metric character name of metric
#' @param value numeric result of applying metric to observed data and fitted distribution
#' @param dens function which takes a vector of coordinates, a list of parameters and
#' a logical determining whether the output should be scaled; returns the fitted distribution density
#' @param seg_start integer histogram index start of uniform segment
#' @param seg_end integer histogram index end of uniform segment
#' @param ... additional parameters, allowing child classes to be built
#'
#' @return an HZModelFit object
new_HZModelFit <- function(
    par = NULL,
    dist = NULL,
    metric = NULL,
    value = NULL,
    dens = NULL,
    seg_start = NULL,
    seg_end = NULL,
    ...
  ){

  # Checking types
  stopifnot(is.list(par) || is.null(par))
  stopifnot(is.character(dist))
  stopifnot(is.character(metric))
  stopifnot(is.numeric(value) || is.na(value))
  stopifnot(is.function(dens))
  stopifnot(is.integer(seg_start) || is.null(seg_start))
  stopifnot(is.integer(seg_end) || is.null(seg_start))

  # Creating object
  x <- list(
    par = par,
    dist = dist,
    metric = metric,
    value = value,
    dens = dens,
    seg_start = seg_start,
    seg_end = seg_end,
    ...
  )
  class(x) <- "HZModelFit"

  return(x)
}

#' @export
print.HZModelFit <- function(x, ...){

  # Distribution/Metric/Observations
  cat("ModelFit\n")
  cat("Distribution: ", x$dist, "\n")
  if(x$metric == "mle"){
    cat("MLE: ", x$value, "\n\n")
  } else {
    cat("Metric: ", x$metric, "\n")
    cat("Metric(obs, exp): ", x$value, "\n\n")
  }

  # Parameters
  if(!is.null(x$par)){
    cat("Fitted parameters:\n")

    param_names <- names(x$par)
    param_vals <- unlist(x$par)
    param_vals <- as.character(formatC(param_vals, digits = 4))
    spacing <- max(nchar(c(param_names, param_vals)))

    param_names <- stringr::str_pad(param_names, width = spacing)
    param_vals <- stringr::str_pad(param_vals, width = spacing)

    cat(param_names, "\n")
    cat(param_vals, "\n")
  }

  # Truncated Uniform
  if(!is.null(x$seg_start) && !is.null(x$seg_end)){
    cat("Fitted segment:\n")

    start_end <- as.character(c(x$seg_start, x$seg_end))
    spacing <- max(nchar(c("Start", "End", start_end)))
    cat(stringr::str_pad(c("Start", "End"), width = spacing), "\n")
    cat(stringr::str_pad(start_end, width = spacing), "\n")
  }

  # Extra line
  cat("\n")

  return(invisible(x))
}

# TODO:
# A function that retrieves density function?
# A function that retrieves other attributes?
