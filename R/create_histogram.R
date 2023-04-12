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

#' @exportS3Method create_histogram Histogram
create_histogram.Histogram <- function(
    x,
    type = c('percent', 'density', 'count'),
    col = 'white',
    xat = unique(c(x$interval_start, x$interval_end)),
    yat = NA,
    xlab.label = '',
    ylab.label = type,
    ylimits = c(0, max(plot_data$y) + 0.05 * max(plot_data$y)),
    xlimits = range(c(x$interval_start, x$interval_end)) + c(-1, 1)
    ) {
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
    bin_width = x$interval_end - x$interval_start,
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
