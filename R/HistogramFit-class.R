
#' @export
print.HistogramFit <- function(x, ...){

  print.Histogram(x) # NextMethod?

  # Print results
  cat("\nResults\n")
  cat("\tNumber of segments: ", length(x$models), "\n")

  # Print parameters
  cat("\nParameters\n")
  cat("\teps: ", x$eps, "\n")
  cat("\tmetrics: ", paste0(x$histogram_metric, collapse = ", "), "\n")
  cat("\tdistributions: ", paste0(x$distributions, collapse = ", "), "\n")
  cat("\tremove low entropy: ", x$remove_low_entropy, "\n")
  cat("\tmaximumize uniform fitting: ", x$max_uniform, "\n")
  cat("\ttruncated models: ", x$truncated_models, "\n")

  invisible(x)
}
