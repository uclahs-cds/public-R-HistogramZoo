
#' @export
print.HistogramFit = function(x, ...){

  print.Histogram(x) # NextMethod?

  # Print results
  cat("\nResults\n")
  cat("\tNumber of segments: ", length(x$models), "\n")

  # Print parameters
  cat("\nParameters\n")
  cat("\teps: ", x$eps, "\n")
  cat("\tmetrics: ", paste0(x$histogram.metric, collapse = ", "), "\n")
  cat("\tdistributions: ", paste0(x$distributions, collapse = ", "), "\n")
  cat("\tremove low entropy: ", x$remove.low.entropy, "\n")
  cat("\tmaximumize uniform fitting: ", x$max.uniform, "\n")
  cat("\ttruncated models: ", x$truncated.models, "\n")

  invisible(x)
}
