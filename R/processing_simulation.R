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

#' Intervals [a1, b1], [a2,b2]
#' Need a <= b, c <= d
#' @return TRUE if [a1,b1] and [a2,b2] overlap
int_overlaps <- function(a1, b1, a2, b2) {
  stopifnot(a1 <= a2 & b1 <= b2)
  b1 <= a2 & a1 <= b2
}

#' Not space efficient at all but faster (?) for integers
intervals_overlap_union <- function(intervals1, intervals2) {
  int1_merged <- unique(unlist(lapply(intervals1, function(x) seq(from = x[1], to = x[2]))))
  int2_merged <- unique(unlist(lapply(intervals2, function(x) seq(from = x[1], to = x[2]))))
  int_intersect <- length(intersect(int1_merged, int2_merged))
  int_union <- length(unique(c(int1_merged, int2_merged)))
  list(
    intersect = int_intersect,
    union = int_union,
    jaccard = int_intersect / int_union
  )
}

#' Computes jaccard for `overlap` data from multi-peak simulation
compute_overall_jaccard <- function(ex) {
  actual_ex <- ex[!duplicated(ex$peak_num), ]
  actual_peaks <- mapply(c, actual_ex$peak_min, actual_ex$peak_max, SIMPLIFY = F)

  fit_ex <- ex[!duplicated(ex$fit_peak_num), ]
  fit_peaks <- mapply(c, fit_ex$fit_peak_start, fit_ex$fit_peak_end, SIMPLIFY = F)

  intervals_overlap_union(actual_peaks, fit_peaks)
}
