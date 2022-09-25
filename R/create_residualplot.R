

#' create_residualplot creates a residual plot between fitted and observed data
#'
#' @param histogram_obj a `Histogram` or `HistogramFit` object
#' @param model_name one of the metrics used to fit models (e.g. Jaccard) and "consensus", default consensus
#' @param add_changepoint_lines logical, whether vertical lines should be drawn where residuals change from positive to negative and vice versa, default FALSE
#' @param abline.h allow horizontal line to be drawn, default to 0
#' @param abline.v allow vertical line to be drawn, default to NULL
#' @param abline.col line colour, defaults to lightgrey
#' @param abline.lwd line width, defaults to 0.01
#' @param abline.lty line style, defaults to dotted
#' @inheritParams BoutrosLab.plotting.general::create.scatterplot
#'
#' @return residual scatterplot, a Trellis object. For further details, see the 'Lattice' R package.
#' @export
#'
#' @examples \dontrun{
#' x = Histogram(c(0, 0, 1, 2, 3, 2, 1, 2, 3, 4, 5, 3, 1, 0))
#' results = segment_and_fit(x, eps = 0.005)
#' create_residualplot(results)
#' }
create_residualplot <- function(
    histogram_obj, model_name,
    add_changepoint_lines,
    abline.h, abline.v, abline.col, abline.lwd, abline.lty,
    filename,
    main, main.just, main.x, main.y, main.cex,
    xlab.label, ylab.label, xlab.cex, ylab.cex, xlab.col, ylab.col,
    xlab.top.label, xlab.top.cex, xlab.top.col, xlab.top.just, xlab.top.x, xlab.top.y,
    xlimits = c(histogram_obj$interval_start[1], histogram_obj$interval_end[length(histogram_obj)]), 
    ylimits, xat, yat, xaxis.lab, yaxis.lab, xaxis.cex, yaxis.cex,
    xaxis.rot, yaxis.rot, xaxis.fontface, yaxis.fontface, xaxis.col, yaxis.col, xaxis.tck, yaxis.tck,
    type, cex, pch, col, col.border, lty, lwd, alpha,
    axes.lwd,
    key, legend,
    top.padding, bottom.padding, right.padding, left.padding,
    key.top, key.left.padding, ylab.axis.padding, axis.key.padding,
    x.spacing, y.spacing,
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
  model_name = c("consensus", histogram_obj$metric),
  add_changepoint_lines = F,
  abline.h = 0,
  abline.v = NULL,
  abline.col = "lightgrey",
  abline.lwd = 0.01,
  abline.lty = "dotted",
  filename = NULL,
  main = histogram_obj$region_id, main.just = 'center', main.x = 0.5, main.y = 0.5, main.cex = 3,
  xlab.label = if(inherits(histogram_obj, "GenomicHistogram")) histogram_obj$chr else "Interval Coordinates",
  ylab.label = "Observed - Fitted",
  xlab.cex = 2, ylab.cex = 2,xlab.col = 'black', ylab.col = 'black', xlab.top.label = NULL, xlab.top.cex = 2, xlab.top.col = 'black',
  xlab.top.just = 'center', xlab.top.x = 0.5, xlab.top.y = 0,
  xlimits = c(histogram_obj$interval_start[1], histogram_obj$interval_end[length(histogram_obj)]),
  ylimits = NULL, xat = TRUE, yat = TRUE, xaxis.lab = NA, yaxis.lab = NA, xaxis.cex = 1.5, yaxis.cex = 1.5,
  xaxis.rot = 0, yaxis.rot = 0, xaxis.fontface = 'bold', yaxis.fontface = 'bold', xaxis.col = 'black', yaxis.col = 'black', xaxis.tck = c(1,1), yaxis.tck = c(1,1),
  type = 'p', cex = 0.75, pch = 19, col = 'black', col.border = 'black', lty = 1, lwd = 1, alpha = 1,
  axes.lwd = 1,
  key = list(text = list(lab = c(''))),
  legend = NULL,
  top.padding = 0.1, bottom.padding = 0.7, right.padding = 0.1, left.padding = 0.5,
  key.top = 0.1, key.left.padding = 0, ylab.axis.padding = 1, axis.key.padding = 1,
  x.spacing = 0, y.spacing = 0,
  add.rectangle = FALSE, xleft.rectangle = NULL, ybottom.rectangle = NULL, xright.rectangle = NULL, ytop.rectangle = NULL, col.rectangle = 'transparent', alpha.rectangle = 1,
  add.points = FALSE, points.x = NULL, points.y = NULL, points.pch = 19, points.col = 'black', points.col.border = 'black', points.cex = 1,
  add.line.segments = FALSE, line.start = NULL, line.end = NULL, line.col = 'black', line.lwd = 1,
  add.text = FALSE, text.labels = NULL, text.x = NULL, text.y = NULL, text.col = 'black', text.cex = 1, text.fontface = 'bold',
  height = 6, width = 6, size.units = 'in', resolution = 1600, enable.warnings = FALSE,
  ...
){

  # Error checking
  stopifnot(inherits(histogram_obj, "HistogramFit"))
  model_name <- match.arg(model_name, c("consensus", histogram_obj$metric))
  stopifnot(is.logical(add_changepoint_lines))

  # Extracting histogram_data
  histogram_data <- histogram_obj$histogram_data
  # choosing the midpoint of the start/end as the label
  labels_x <- rowMeans(cbind(histogram_obj$interval_start, histogram_obj$interval_end))
  plotting_data <- data.frame("density" = histogram_data, "labels_x" = labels_x)

  # Distribution fit data
  mods <- lapply(histogram_obj$models, `[[`, model_name)
  distribution_plotting_data <- lapply(mods, function(m) {
    x <- seq(m$seg_start, m$seg_end, by = 1)
    dens <- m$dens(x = seq_along(x), mpar = m$par)
    return(
      data.frame("fitted" = dens, "labels_x" = labels_x[x])
    )
  })
  distribution_plotting_data <- do.call('rbind.data.frame', distribution_plotting_data)

  # Calculating residuals
  plotting_data <- merge(plotting_data, distribution_plotting_data, by = "labels_x", all = T)
  plotting_data$Residuals <- plotting_data$density - plotting_data$fitted

  # Add changepoint lines
  if(add_changepoint_lines){
    abline.v <- unique(
      c(abline.v, which(abs(diff(sign(plotting_data$Residuals))) == 2))
    )
    abline.h <- unique(
      c(abline.h, 0)
    )
  }

  # Plotting
  plt <-  BoutrosLab.plotting.general::create.scatterplot(
    Residuals ~ labels_x,
    plotting_data,
    filename = filename,
    # Lines
    abline.h = abline.h,
    abline.v = abline.v,
    abline.lty = abline.lty,
    abline.col = abline.col,
    abline.lwd = abline.lwd,
    # Titles and labels
    main = main,
    main.just = main.just,
    main.x = main.x,
    main.y = main.y,
    main.cex = main.cex,
    xlab.label = xlab.label,
    ylab.label = ylab.label,
    xlab.cex = xlab.cex,
    ylab.cex = ylab.cex,
    xlab.col = xlab.col,
    ylab.col = ylab.col,
    xlab.top.label = xlab.top.label,
    xlab.top.cex = xlab.top.cex,
    xlab.top.col = xlab.top.col,
    xlab.top.just = xlab.top.just,
    xlab.top.x = xlab.top.x,
    xlab.top.y = xlab.top.y,
    xlimits = xlimits,
    ylimits = ylimits,
    xat = xat,
    yat = yat,
    xaxis.lab = xaxis.lab,
    yaxis.lab = yaxis.lab,
    xaxis.cex = xaxis.cex,
    yaxis.cex = yaxis.cex,
    xaxis.rot = xaxis.rot,
    yaxis.rot = yaxis.rot,
    xaxis.fontface = xaxis.fontface,
    yaxis.fontface = yaxis.fontface,
    xaxis.col = xaxis.col,
    yaxis.col = yaxis.col,
    xaxis.tck = xaxis.tck,
    yaxis.tck = yaxis.tck,
    type = type,
    cex = cex,
    pch = pch,
    col = col,
    col.border = col.border,
    lty = lty,
    lwd = lwd,
    alpha = alpha,
    axes.lwd = axes.lwd,
    key = key,
    legend = legend,
    top.padding = top.padding,
    bottom.padding = bottom.padding,
    right.padding = right.padding,
    left.padding = left.padding,
    key.top = key.top,
    key.left.padding = key.left.padding,
    ylab.axis.padding = ylab.axis.padding,
    axis.key.padding = axis.key.padding,
    x.spacing = x.spacing,
    y.spacing = y.spacing,
    add.rectangle = add.rectangle,
    xleft.rectangle = xleft.rectangle,
    ybottom.rectangle = ybottom.rectangle,
    xright.rectangle = xright.rectangle,
    ytop.rectangle = ytop.rectangle,
    col.rectangle = col.rectangle,
    alpha.rectangle = alpha.rectangle,
    add.points = add.points,
    points.x = points.x,
    points.y = points.y,
    points.pch = points.pch,
    points.col = points.col,
    points.col.border = points.col.border,
    points.cex = points.cex,
    add.line.segments = add.line.segments,
    line.start = line.start,
    line.end = line.end,
    line.col = line.col,
    line.lwd = line.lwd,
    add.text = add.text,
    text.labels = text.labels,
    text.x = text.x,
    text.y = text.y,
    text.col = text.col,
    text.cex = text.cex,
    text.fontface = text.fontface,
    height = height,
    width = width,
    size.units = size.units,
    resolution = resolution,
    enable.warnings = enable.warnings,
    ...
  )

  # Return plot
  return(plt)
}
