
# Defining plotting parameters
distributions <- c(
  "norm",
  "gamma",
  "gamma_flip",
  "unif"
)

distribution_colours <- c(
  "norm" = "darkorange",
  "gamma" = "chartreuse4",
  "gamma_flip" = "chartreuse3",
  "unif" = "darkorchid4"
)

distribution_names <- c(
  "norm" = "Normal",
  "gamma" = "Gamma",
  "gamma_flip" = "Gamma Flip",
  "unif" = "Uniform"
)

#' create_coverageplot creates a histogram/coverage plot with the potential to add annotations
#'
#' @param histogram_obj a `Histogram` or `HistogramFit` object
#' @param model_name if `histogram_obj` is a `HistogramFit` object, one of the metrics used to fit models and "consensus" will be used to plot distribution fit data, default consensus
#' @param add.points allow additional points to be drawn, if `histogram_obj` is a `HistogramFit` object, default TRUE
#' @param points.x The x co-ordinates where additional points should be drawn, if `histogram_obj` is a `HistogramFit` object, default segment_and_fit points
#' @param points.y The y co-ordinates where additional points should be drawn, if `histogram_obj` is a `HistogramFit` object, default segment_and_fit points
#' @param col colour for coverage, default `black`
#' @param lwd lwd for coverage, default 1
#' @param lty lwd for coverage, default 1
#' @param col_distributions colour vector for the distributions, colours should be indicated in the order of `histogram_obj$distributions`, default see `HistogramZoo:::distribution_colours`
#' @param lwd_distributions lwd vector for the distributions, default 2.5
#' @param lty_distributions lty vector for the distributions, default 1
#' @param type `type` in R graphics. Default: if histogram length < 50, plot `h`, otherwise `l`
#' @inheritParams BoutrosLab.plotting.general::create.scatterplot
#'
#' @return Coverage plot, a Trellis object. For further details, see the 'Lattice' R package.
#' @export
#'
#' @examples \dontrun{
#' x = rnorm(10000, mean = 100, sd = 50)
#' x = observations_to_histogram(round(x), histogram_bin_width = 5)
#' results = segment_and_fit(x, eps = 1)
#' create_coverageplot(
#'   results
#' )
#' }
create_coverageplot <- function(
    histogram_obj, model_name,
    col, lwd, lty,
    col_distributions, lwd_distributions, lty_distributions,
    type,
    add.points, points.x, points.y, points.pch, points.col, points.col.border, points.cex,
    filename,
    main, main.just, main.x, main.y, main.cex,
    xlab.label, ylab.label, xlab.cex, ylab.cex, xlab.col, ylab.col,
    xlab.top.label, xlab.top.cex, xlab.top.col, xlab.top.just, xlab.top.x, xlab.top.y,
    xlimits, ylimits, xat, yat, xaxis.lab, yaxis.lab, xaxis.cex, yaxis.cex,
    xaxis.rot, yaxis.rot, xaxis.fontface, yaxis.fontface, xaxis.col, yaxis.col, xaxis.tck, yaxis.tck,
    cex, col.border, pch, alpha,
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
    histogram_obj, model_name = NULL,
    col = 'black', lwd = 1, lty = 1,
    col_distributions = NULL, lwd_distributions = NULL, lty_distributions = NULL,
    type = if(length(histogram_obj) < 50) 'h' else 'l',
    add.points = FALSE, points.x = NULL, points.y = NULL, points.pch = 19, points.col = 'black', points.col.border = 'black', points.cex = 1,
    filename = NULL,
    main = histogram_obj$region_id, main.just = 'center', main.x = 0.5, main.y = 0.5, main.cex = 3,
    xlab.label = if(inherits(histogram_obj, "GenomicHistogram")) histogram_obj$chr else "Interval Coordinates",
    ylab.label = "Histogram Data",
    xlab.cex = 2, ylab.cex = 2,xlab.col = 'black', ylab.col = 'black', xlab.top.label = NULL, xlab.top.cex = 2, xlab.top.col = 'black',
    xlab.top.just = 'center', xlab.top.x = 0.5, xlab.top.y = 0,
    xlimits = c(1, length(histogram_obj)),
    ylimits = NULL, 
    xat = generate_xlabels(histogram_obj, return_xat = T), 
    yat = TRUE, 
    xaxis.lab = generate_xlabels(histogram_obj), 
    yaxis.lab = NA, xaxis.cex = 1.5, yaxis.cex = 1.5,
    xaxis.rot = 0, yaxis.rot = 0, xaxis.fontface = 'bold', yaxis.fontface = 'bold', xaxis.col = 'black', yaxis.col = 'black', xaxis.tck = c(1,1), yaxis.tck = c(1,1),
    cex = 0.75, col.border = 'black', pch = 19, alpha = 1,
    axes.lwd = 1,
    key = list(text = list(lab = c(''))), legend = NULL,
    top.padding = 0.1, bottom.padding = 0.7, right.padding = 0.1, left.padding = 0.5,
    key.top = 0.1, key.left.padding = 0, ylab.axis.padding = 1, axis.key.padding = 1,
    x.spacing = 0, y.spacing = 0,
    abline.h = NULL, abline.v = NULL, abline.col = 'black', abline.lwd = 1, abline.lty = 1,
    add.rectangle = FALSE, xleft.rectangle = NULL, ybottom.rectangle = NULL, xright.rectangle = NULL, ytop.rectangle = NULL, col.rectangle = 'transparent', alpha.rectangle = 1,
    add.line.segments = FALSE, line.start = NULL, line.end = NULL, line.col = 'black', line.lwd = 1,
    add.text = FALSE, text.labels = NULL, text.x = NULL, text.y = NULL, text.col = 'black', text.cex = 1, text.fontface = 'bold',
    height = 6, width = 6, size.units = 'in', resolution = 1600, enable.warnings = FALSE,
    ...
){

  # Error checking
  stopifnot(inherits(histogram_obj, "Histogram"))

  # Extracting histogram_data
  histogram_data <- histogram_obj$histogram_data
  # choosing the midpoint of the start/end as the label
  plotting_data <- data.frame("dens" = histogram_data, "labels_x" = seq(1, length(histogram_data), 1))

  # Plotting
  plt <- BoutrosLab.plotting.general::create.scatterplot(
    dens ~ labels_x,
    data = plotting_data,
    filename = filename,
    # Lines & PCH
    cex = cex,
    col.border = col.border,
    pch = pch,
    lty = lty,
    alpha = alpha,
    col = col,
    lwd = lwd,
    type = type,
    axes.lwd = axes.lwd,
    # Keys and legends and padding
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
    # lines
    abline.h = abline.h,
    abline.v = abline.v,
    abline.col = abline.col,
    abline.lwd = abline.lwd,
    abline.lty = abline.lty,
    add.rectangle = add.rectangle,
    xleft.rectangle = xleft.rectangle,
    ybottom.rectangle = ybottom.rectangle,
    xright.rectangle = xright.rectangle,
    ytop.rectangle = ytop.rectangle,
    col.rectangle = col.rectangle,
    alpha.rectangle = alpha.rectangle,
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
    # Labels
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
    # Extra plotting parameters
    add.points = add.points,
    points.x = points.x,
    points.y = points.y,
    points.pch = points.pch,
    points.col = points.col,
    points.col.border = points.col.border,
    points.cex = points.cex,
    size.units = size.units,
    resolution = resolution,
    enable.warnings = enable.warnings,
    ...
  )

  # Return plot
  return(plt)
}


#' Returns the segment_and_fit x coordinates of the identified points
#' @param histogram_obj a HistogramFit object
#' @return numeric vector, representing segmentation points
return_x_points <- function(histogram_obj){
  return(
      c(histogram_obj$p[,'start'], histogram_obj$p[,'end'])
  )
}

#' Returns the segment_and_fit y coordinates of the identified points
#' @param histogram_obj a HistogramFit object
#' @return numeric vector, representing histogram values at segmentation points
return_y_points <- function(histogram_obj){
  return(
    histogram_obj$histogram_data[
      c(histogram_obj$p[,'start'], histogram_obj$p[,'end'])
    ]
  )
}

#' @export
create_coverageplot.HistogramFit <- function(
  histogram_obj,
  model_name = c("consensus", histogram_obj$metric),
  col = 'black',
  lwd = 1,
  lty = 1,
  col_distributions = distribution_colours[histogram_obj$distributions],
  lwd_distributions = 2.5,
  lty_distributions = 1,
  type = if(length(histogram_obj) < 50) 'h' else 'l',
  add.points = T,
  points.x = return_x_points(histogram_obj),
  points.y = return_y_points(histogram_obj),
  points.pch = 19,
  points.col = 'red',
  points.col.border = 'black',
  points.cex = 1,
  filename = NULL,
  main = histogram_obj$region_id, main.just = 'center', main.x = 0.5, main.y = 0.5, main.cex = 3,
  xlab.label = if(inherits(histogram_obj, "GenomicHistogram")) histogram_obj$chr else "Interval Coordinates",
  ylab.label = "Histogram Data",
  xlab.cex = 2, ylab.cex = 2,xlab.col = 'black', ylab.col = 'black', xlab.top.label = NULL, xlab.top.cex = 2, xlab.top.col = 'black',
  xlab.top.just = 'center', xlab.top.x = 0.5, xlab.top.y = 0,
  xlimits = c(1, length(histogram_obj)),  
  ylimits = NULL,
  xat = generate_xlabels(histogram_obj, return_xat = T), 
  yat = TRUE, 
  xaxis.lab = generate_xlabels(histogram_obj),
  yaxis.lab = NA, xaxis.cex = 1.5, yaxis.cex = 1.5,
  xaxis.rot = 0, yaxis.rot = 0, xaxis.fontface = 'bold', yaxis.fontface = 'bold', xaxis.col = 'black', yaxis.col = 'black', xaxis.tck = c(1,1), yaxis.tck = c(1,1),
  cex = 0.75, col.border = 'black', pch = 19, alpha = 1,
  axes.lwd = 1,
  key = list(text = list(lab = c(''))),
  legend = list(
    right = list(
      fun = lattice::draw.key,
      args = list(
        key = list(
          points = list(
            col = "black",
            pch = 22,
            cex = 2,
            fill = c(col, col_distributions)
          ),
          text = list(
            lab = c("Coverage", distribution_names[histogram_obj$distributions])
          ),
          padding.text = 3,
          cex = 1
        )
      )
    )
  ),
  top.padding = 0.1, bottom.padding = 0.7, right.padding = 0.1, left.padding = 0.5,
  key.top = 0.1, key.left.padding = 0, ylab.axis.padding = 1, axis.key.padding = 1,
  x.spacing = 0, y.spacing = 0,
  abline.h = NULL, abline.v = NULL, abline.col = 'black', abline.lwd = 1, abline.lty = 1,
  add.rectangle = FALSE, xleft.rectangle = NULL, ybottom.rectangle = NULL, xright.rectangle = NULL, ytop.rectangle = NULL, col.rectangle = 'transparent', alpha.rectangle = 1,
  add.line.segments = FALSE, line.start = NULL, line.end = NULL, line.col = 'black', line.lwd = 1,
  add.text = FALSE, text.labels = NULL, text.x = NULL, text.y = NULL, text.col = 'black', text.cex = 1, text.fontface = 'bold',
  height = 6, width = 6, size.units = 'in', resolution = 1600, enable.warnings = FALSE,
  ...
){

  # Error checking
  stopifnot(inherits(histogram_obj, "HistogramFit"))
  if(length(histogram_obj$distributions) != length(col_distributions)){
    warning("Number of distributions fit does not equal the number of elements in col_distributions.")
  }
  model_name <- match.arg(model_name, c("consensus", histogram_obj$metric))

  # Extracting histogram_data
  histogram_data <- histogram_obj$histogram_data
  # choosing the midpoint of the start/end as the label
  plotting_data <- data.frame(
    "dens" = histogram_data, 
    "labels_x" = seq(1, length(histogram_data), 1), 
    "dist" = "coverage"
  )
  segment_dists <- unlist(lapply(histogram_obj$models, function(x) x$consensus$dist))
  dist_colors <- col_distributions[segment_dists]
  # Distribution fit data
  mods <- lapply(histogram_obj$models, `[[`,  model_name)
  distribution_plotting_data <- lapply(seq_along(mods), function(i) {
    m <- mods[[i]]
    x <- seq(m$seg_start, m$seg_end, by = 1)
    dens <- m$dens(x = seq_along(x), mpar = m$par)
    return(
      data.frame("dens" = dens, "labels_x" = x, "dist" = m$dist, "segment" = i)
    )
  })
  distribution_plotting_data <- do.call('rbind.data.frame', distribution_plotting_data)

  # Factoring plotting data distribution
  distribution_plotting_data$dist <- factor(distribution_plotting_data$dist, levels = histogram_obj$distributions)

  # Plotting
  base_plot <- BoutrosLab.plotting.general::create.scatterplot(
    dens ~ labels_x,
    data = plotting_data,
    filename = filename,
    # Lines & PCH
    cex = cex,
    col.border = col.border,
    pch = pch,
    lty = lty,
    alpha = alpha,
    col = col,
    lwd = lwd,
    type = type,
    axes.lwd = axes.lwd,
    # Keys and legends and padding
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
    # lines
    abline.h = abline.h,
    abline.v = abline.v,
    abline.col = abline.col,
    abline.lwd = abline.lwd,
    abline.lty = abline.lty,
    add.rectangle = add.rectangle,
    xleft.rectangle = xleft.rectangle,
    ybottom.rectangle = ybottom.rectangle,
    xright.rectangle = xright.rectangle,
    ytop.rectangle = ytop.rectangle,
    col.rectangle = col.rectangle,
    alpha.rectangle = alpha.rectangle,
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
    # Labels
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
    # Extra plotting parameters
    add.points = add.points,
    points.x = points.x,
    points.y = points.y,
    points.pch = points.pch,
    points.col = points.col,
    points.col.border = points.col.border,
    points.cex = points.cex,
    size.units = size.units,
    resolution = resolution,
    enable.warnings = enable.warnings,
    ...
  )


  dist_plot <- BoutrosLab.plotting.general::create.scatterplot(
    formula = dens ~ labels_x,
    data = distribution_plotting_data,
    # Groups
    groups = distribution_plotting_data$segment,
    col = dist_colors,
    lwd = lwd_distributions,
    lty = lty_distributions,
    # Labels
    main = main,
    xlab.label = xlab.label,
    ylab.label = ylab.label,
    # Lines & PCH
    type = 'l'
  )


  # Return plot
  return(base_plot + dist_plot)
}
