#' @export
#' @rdname create_histogram
create_histogram <- function(x, ...) {
  UseMethod('create_histogram')
}

#' @exportS3Method density Histogram
density.Histogram <- function(x) {
  counts <- x$histogram_data
  bin_widths <- x$interval_end - x$interval_start
  n <- sum(x$histogram_data)
  counts / (n * bin_widths)
}

#' Creates a histogram with potentially variable length bins
#'
#' @param x Histogram
#'
#' @param type 'percent', 'density' or 'count'. Note that count or percentage are misleading if unequal bin-widths
#' @param col Fill color of the histogram
#' @param xat xat
#' @param yat yat
#' @param xlab.label x-axis label
#' @param ylab.label y-axis label
#' @param ylimits ylim in lattice
#' @param xlimits xlim in lattice
#' @param ... additional arguments to be pass into lattice::xyplot
#'
#' @exportS3Method create_histogram Histogram
#' @rdname create_histogram
#' @examples
#' # Equal length bin widths
#' my_histogram <- Histogram(
#'  histogram_data = c(1,2,3,1,5,6),
#'  interval_start = c(0,1,2,3,4,5),
#'  interval_end = c(1,2,3,4,5,6),
#'  bin_width = 1,
#'  region_id = "my_histogram"
#' )
#'
#' create_histogram(my_histogram, type = 'count')
#'
#' my_histogram$interval_end[6] <- 10
#'
#' # Defaults to density if unequal bin-widths
#' create_histogram(my_histogram)
create_histogram.Histogram <- function(
    x,
    type = c('percent', 'density', 'count'),
    col = 'white',
    xat = unique(c(x$interval_start, x$interval_end)),
    yat = NA,
    xlab.label = '',
    ylab.label = type,
    ylimits = c(0, max(plot_data$y) + 0.05 * max(plot_data$y)),
    xlimits = range(c(x$interval_start, x$interval_end)) + c(-1, 1),
    ...
    ) {
  bin_width <- x$interval_end - x$interval_start
  if (length(unique(bin_width)) > 1) {
    # Default to density if unequal bin-widths
    if (missing(type)) {
      type <- 'density'
      } else if (type %in% c('percent', 'count')) {
      warning('Using `percentage` or `count` with unequal bin-widths may be misleading as AREAS are not proportional to probability.')
      }
    }
  type <- match.arg(type)
  plot_data <- list(
    y = switch(
      type,
      count = x$histogram_data,
      percent = 100 * x$histogram_data/sum(x$histogram_data),
      density = density(x)
      ),
    x = x$interval_start
  )

  lattice::xyplot(
    y ~ x,
    data = plot_data,
    bin_width = bin_width,
    ylim = ylimits,
    xlim = xlimits,
    scales = list(
      x = list(
        at = xat,
        fontface = 'bold'
        ),
      y = list(
        fontface = 'bold',
        at = yat
      )
    ),
    xlab = list(
      label = xlab.label,
      fontface = 'bold'
      ),
    ylab = list(
      label = ylab.label,
      fontface = 'bold'
      ),
    ...,
    panel = function(x, y, bin_width = 1) {
        lattice::panel.rect(
            x = x,
            y = 0,
            height = y,
            width = bin_width,
            just = c("left", "bottom"),
            col = col
          )
    })
}


create_histogram.HistogramFit <- function(
    x,
    ...) {
  create_histogram.Histogram(
    x,
    ...
    )
}
