#' A helper function that generates labels at specific Histogram bins
#'
#' @param histogram_obj a `Histogram` or `HistogramFit` object
#' @param xat the indices of the Histogram at which to generate labels
#'
#' @return labels at `xat` coordinates
labels_helper <- function(
    histogram_obj,
    xat
){

  interval_start_type <- c('[', rep('(', length(histogram_obj$interval_start) - 1))
  interval_start <- paste0(interval_start_type, histogram_obj$interval_start)
  interval_end <- paste0(histogram_obj$interval_end, ']')

  labels <- ifelse(
    histogram_obj$interval_start == histogram_obj$interval_end,
    as.character(histogram_obj$interval_start),
    paste0(interval_start, ",", interval_end)
  )

  return(labels[xat])
}

#' Generates x labels for histogram-based plots, `create_coverageplot` and `create_residualplot`
#'
#' @param histogram_obj a `Histogram` or `HistogramFit` object
#' @param n_labels number of x axis labels
#' @param return_xat logical, whether to return the indices of the labels, default FALSE, returns the labels
#' @param disjoint logical, whether to check if the `histogram_obj` contains disjoint intervals. If found only breakpoints will be labelled
#' if no disjoint intervals are found, x labels default to standard
#'
#' @return a vector of labels (character) or indices (numeric)
#' @export
#'
#' @examples \dontrun{
#'  x <- observations_to_histogram(rnorm(10000))
#'  generate_xlabels(x, n_labels = 5, xat = F, disjoint = F)
#' }
generate_xlabels <- function(
    histogram_obj,
    n_labels = 5,
    return_xat = F,
    disjoint = T
){

  xat = c()
  if(disjoint){
    bin_width <- histogram_obj$bin_width
    xat <- which(diff(histogram_obj$interval_start) > bin_width) # This chooses the left coordinate
    if(length(xat) > 0){
      x_labels <- paste0(
        labels_helper(histogram_obj, xat), "|",
        labels_helper(histogram_obj, xat + 1)
      )
      if(!1 %in% xat){
        xat <- c(1, xat)
        x_labels <- c(labels_helper(histogram_obj, 1), x_labels)
      }
      if(!length(histogram_obj) %in% (xat + 1)){
        xat <- c(xat, length(histogram_obj))
        x_labels <- c(x_labels, labels_helper(histogram_obj, length(histogram_obj)))
      }
    }
  }
  if(length(xat) == 0){
    xat <- seq(from = 1, to = length(histogram_obj), length.out = n_labels)
    xat <- round(xat)
    x_labels <- labels_helper(histogram_obj, xat)
  }

  if(return_xat){
    return(xat)
  } else {
    return(x_labels)
  }

}
