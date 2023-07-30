# Should be un-scale
scale.param <- function(x, param.range) {
  (param.range[2] - param.range[1]) * x + param.range[1]
}


scale.param2 <- function(x, param.range) {
 (x -  param.range[1]) / (param.range[2] - param.range[1])
}

summarize_results_error <- function(x) {
  if (! is(x,  'try-error')) {
        models <- x$models
        metrics <- x$metric
        all_metric_results <- do.call(plyr::rbind.fill, lapply(metrics, summarize_results, result = x))

        rtn <- all_metric_results
        rtn$num_segments <- length(models)
        rtn
    } else {
        rtn <- data.frame(error = x[[1]])
    }
  return(rtn);
}

#' @export
overlap_size <- function(a1, a2, b1, b2) {
  max(0, min(a2, b2) - max(a1, b1))
}

#' @export
union_size <- function(a1, a2, b1, b2) {
  if (overlap_size(a1, a2, b1, b2) > 0) max(a2, b2) - min(a1, b1)
  else (a2 - a1) + (b2 - b1)
}
