
# Defining plotting parameters
distributions = c(
  "coverage",
  "norm",
  "gamma",
  "unif"
)

distribution_colours = c(
  "coverage" = "black",
  "norm" = "darkorange",
  "gamma" = "chartreuse4",
  "unif" = "darkorchid4"
)

distribution_lwd = c(
  "coverage" = 1,
  "norm" = 2.5,
  "gamma" = 2.5,
  "unif" = 2.5
)

#' create_coverageplot creates a histogram/coverage plot with the potential to add annotations
#'
#' @param histogram_obj a `Histogram` or `HistogramFit` object
#' @param model_name if `histogram_obj` is a `HistogramFit` object, one of the metrics used to fit models and "consensus" will be used to plot distribution fit data, default consensus
#' @param add.points allow additional points to be drawn, if `histogram_obj` is a `HistogramFit` object, default TRUE
#' @param points.x The x co-ordinates where additional points should be drawn, if `histogram_obj` is a `HistogramFit` object, default segment_and_fit points
#' @param points.y The y co-ordinates where additional points should be drawn, if `histogram_obj` is a `HistogramFit` object, default segment_and_fit points
#' @param col colour vector for the distributions, default see `HistogramZoo:::distribution_colours`
#' @param lwd lwd vector for the distributions, default see `HistogramZoo:::distribution_lwd`
#' @inheritParams BoutrosLab.plotting.general::create.scatterplot 
#' 
#' @return Coverage plot, a Trellis object. For further details, see the 'Lattice' R package.
#' @export
#'
#' @examples \dontrun{
#' x = Histogram(c(0, 0, 1, 2, 3, 2, 1, 2, 3, 4, 5, 3, 1, 0))
#' results = segment_and_fit(x, eps = 0.005)
#' create_coverageplot(results)
#' }
create_coverageplot <- function(
    histogram_obj, model_name,
    col, lwd,
    add.points, points.x, points.y, points.pch, points.col, points.col.border, points.cex,
    filename,
    main, main.just, main.x, main.y, main.cex,
    xlab.label, ylab.label, xlab.cex, ylab.cex, xlab.col, ylab.col,
    xlab.top.label, xlab.top.cex, xlab.top.col, xlab.top.just, xlab.top.x, xlab.top.y,
    xlimits, ylimits, xat, yat, xaxis.lab, yaxis.lab, xaxis.cex, yaxis.cex,
    xaxis.rot, yaxis.rot, xaxis.fontface, yaxis.fontface, xaxis.col, yaxis.col, xaxis.tck, yaxis.tck,
    type, cex, pch, lty, alpha,
    axes.lwd, 
    key, legend,
    top.padding, bottom.padding, right.padding, left.padding,
    key.top, key.left.padding, ylab.axis.padding, axis.key.padding,
    x.spacing, y.spacing,  
    abline.h, abline.v, abline.col, abline.lwd, abline.lty,
    add.rectangle, xleft.rectangle, ybottom.rectangle, xright.rectangle, ytop.rectangle, col.rectangle, alpha.rectangle,
    add.line.segments, line.start, line.end, line.col, line.lwd,
    add.text, text.labels, text.x, text.y, text.col, text.cex, text.fontface,
    height, width, size.units, resolution, enable.warnings, 
    ...
){
  UseMethod('create_coverageplot')
}


#' @export
create_coverageplot.Histogram <- function(
    histogram_obj,
    main = histogram_obj$region_id,
    ylab.label = "Histogram Data",
    xlab.label = if(inherits(histogram_obj, "GenomicHistogram")) histogram_obj$chr else "Interval Coordinates",
    ...
){
  
  # Error checking
  stopifnot(inherits(histogram_obj, "Histogram"))

  # Extracting histogram_data
  x <- histogram_obj$histogram_data
  # choosing the midpoint of the start/end as the label
  labels_x <- rowMeans(cbind(histogram_obj$interval_start, histogram_obj$interval_end))
  plotting.data <- data.frame("x" = x, "labels.x" = labels_x)

  # Plotting
  plt <- BoutrosLab.plotting.general::create.scatterplot(
    x ~ labels.x,
    data = plotting.data,
    # Lines & PCH
    type = c('a'),
    # Labels
    main = main,
    xlab.label = xlab.label,
    ylab.label = ylab.label,
    # Extra plotting parameters
    ...
  )
  
  # Return plot
  return(plt)
}


#' Returns the segment_and_fit x coordinates of the identified points
return_x_points = function(histogram_obj){
  labels_x <- rowMeans(cbind(histogram_obj$interval_start, histogram_obj$interval_end))
  return(
    labels_x[
      c(histogram_obj$p[,'start'], histogram_obj$p[,'end'])
      ]
  )
}

#' Returns the segment_and_fit y coordinates of the identified points
return_y_points = function(histogram_obj){
  return(
    histogram_obj$histogram_data[
      c(histogram_obj$p[,'start'], histogram_obj$p[,'end'])
    ]
  )
}

#' @export
create_coverageplot.HistogramFit <- function(
  histogram_obj,
  model_name = c("consensus", histogram_obj$histogram_metric),
  col = distribution_colours,
  lwd = distribution_lwd,
  main = histogram_obj$region_id,
  ylab.label = "Histogram Data",
  xlab.label = if(inherits(histogram_obj, "GenomicHistogram")) histogram_obj$chr else "Interval Coordinates",
  add.points = T,
  points.x = return_x_points(histogram_obj),
  points.y = return_y_points(histogram_obj),
  points.pch = 19,
  points.col = 'red',
  ...
){

  # Error checking
  stopifnot(inherits(histogram_obj, "HistogramFit"))
  model_name <- match.arg(model_name, c("consensus", histogram_obj$histogram_metric))

  # Extracting histogram_data
  x <- histogram_obj$histogram_data
  # choosing the midpoint of the start/end as the label
  labels_x <- rowMeans(cbind(histogram_obj$interval_start, histogram_obj$interval_end))
  plotting.data <- data.frame("x" = x, "labels.x" = labels_x, "dist" = "coverage")

  # distribution fit data
  mods <- lapply(histogram_obj$models, `[[`,  model_name)
  distribution_plotting_data <- lapply(mods, function(m) {
    x <- seq(m$seg.start, m$seg.end, by = 1)
    dens <- m$dens(x = seq_along(x), mpar = m$par)
    return(
      data.frame("x" = dens, "labels.x" = labels_x[x], "dist" = m$dist)
    )
  })
  distribution_plotting_data <- do.call('rbind.data.frame', distribution_plotting_data)
  plotting.data <- rbind(plotting.data, distribution_plotting_data)

  # Factoring plotting data distribution
  plotting.data$dist <- factor(plotting.data$dist, levels = distributions)

  # Plotting
  plt <- BoutrosLab.plotting.general::create.scatterplot(
    x ~ labels.x,
    data = plotting.data,
    # Groups
    groups = plotting.data$dist,
    col = col,
    lwd = lwd,
    # Labels
    main = main,
    xlab.label = xlab.label,
    ylab.label = ylab.label,
    # Lines & PCH
    type = c('a'),
    # Adding extra points
    add.points = add.points,
    points.x = points.x,
    points.y = points.y,
    points.pch = points.pch,
    points.col = points.col,
    # Extra plotting parameters
    ...
  )

  # Return plot
  return(plt)
}
