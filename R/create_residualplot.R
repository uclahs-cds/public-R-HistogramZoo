

#' create_residualplot creates a residual plot between fitted and observed data
#'
#' @param histogram_obj a `Histogram` or `HistogramFit` object
#' @param model_name One of the metrics used to fit models (e.g. Jaccard) and "consensus", default consensus
#' @inheritParams BoutrosLab.plotting.general::create.scatterplot 
#'
#' @return Residual scatterplot, a Trellis object. For further details, see the 'Lattice' R package.
#' @export
#'
#' @examples \dontrun{
#' x = Histogram(c(0, 0, 1, 2, 3, 2, 1, 2, 3, 4, 5, 3, 1, 0))
#' results = segment_and_fit(x, eps = 0.005)
#' create_residualplot(results)
#' }
create_residualplot <- function(
    histogram_obj, model_name,
    filename,
    main, main.just, main.x, main.y, main.cex,
    xlab.label, ylab.label, xlab.cex, ylab.cex, xlab.col, ylab.col,
    xlab.top.label, xlab.top.cex, xlab.top.col, xlab.top.just, xlab.top.x, xlab.top.y,
    xlimits, ylimits, xat, yat, xaxis.lab, yaxis.lab, xaxis.cex, yaxis.cex,
    xaxis.rot, yaxis.rot, xaxis.fontface, yaxis.fontface, xaxis.col, yaxis.col, xaxis.tck, yaxis.tck,
    type, cex, pch, col, col.border, lty, lwd, alpha,
    axes.lwd, 
    key, legend,
    top.padding, bottom.padding, right.padding, left.padding,
    key.top, key.left.padding, ylab.axis.padding, axis.key.padding,
    x.spacing, y.spacing,  
    abline.h, abline.v, abline.col, abline.lwd, abline.lty,
    add.rectangle, xleft.rectangle, ybottom.rectangle, xright.rectangle, ytop.rectangle, col.rectangle, alpha.rectangle,
    add.points, points.x, points.y, points.pch, points.col, points.col.border, points.cex,
    add.line.segments, line.start, line.end, line.col, line.lwd,
    add.text, text.labels, text.x, text.y, text.col, text.cex, text.fontface,
    height, width, size.units, resolution, enable.warnings, 
    ...
){
  UseMethod('create_residualplot')
}

#' @export
create_residualplot.HistogramFit = function(
  histogram_obj,
  model_name = c("consensus", histogram_obj$histogram_metric),
  ...
){

  # Error checking
  stopifnot(inherits(histogram_obj, "HistogramFit"))
  model_name <- match.arg(model_name, c("consensus", histogram_obj$histogram_metric))

  # Extracting histogram_data
  x <- histogram_obj$histogram_data
  # choosing the midpoint of the start/end as the label
  labels_x <- rowMeans(cbind(histogram_obj$interval_start, histogram_obj$interval_end))
  plotting.data <- data.frame("density" = x, "labels.x" = labels_x)

  # Distribution fit data
  mods <- lapply(histogram_obj$models, `[[`, model_name)
  distribution_plotting_data <- lapply(mods, function(m) {
    x <- seq(m$seg.start, m$seg.end, by = 1)
    dens <- m$dens(x = seq_along(x), mpar = m$par)
    return(
      data.frame("fitted" = dens, "labels.x" = labels_x[x])
    )
  })
  distribution_plotting_data <- do.call('rbind.data.frame', distribution_plotting_data)

  # Calculating residuals
  plotting.data <- merge(plotting.data, distribution_plotting_data, by = "labels.x", all = T)
  plotting.data$Residuals <- plotting.data$density - plotting.data$fitted

  # Plotting
  plt <-  BoutrosLab.plotting.general::create.scatterplot(
    Residuals ~ labels.x,
    plotting.data,
    ...
  )

  # Return plot
  return(plt)
}
