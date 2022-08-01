
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
#' @param model_name One of the metrics used to fit models (e.g. Jaccard) and "consensus", default consensus
#' @param ... Additional parameters for the BoutrosLab.plotting.general function create.scatterplot
#'
#' @return Coverage plot, a Trellis object. For further details, see the 'Lattice' R package.
#' @export
#'
#' @examples
create_coverageplot <- function(histogram_obj, model_name){
  UseMethod('create_coverageplot')
}

#' @export
create_coverageplot.HistogramFit <- function(
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
  plotting.data <- data.frame(
    "x" = x,
    "labels.x" = 1:length(x),
    "dist" = "coverage")

  # points
  points.x <- labels_x[c(histogram_obj$p[,'start'], histogram_obj$p[,'end'])]
  points.y <- x[points.x]

  # distribution
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
    col = distribution_colours,
    lwd = distribution_lwd,
    # Lines & PCH
    type = c('a'),
    # Adding extra points
    add.points = T,
    points.x = points.x,
    points.y = points.y,
    points.pch = 19,
    points.col = 'red',
    # Extra plotting parameters
    ...
  )

  # Return plot
  return(plt)
}
