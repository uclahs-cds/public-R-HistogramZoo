
# Defining plotting parameters
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

#' create.coverageplot creates a histogram/coverage plot with the potential to add annotations
#'
#' @param x A numeric vector representing coverage quantitation. Vector names (optional) represent the indices of the coverage quantitation.
#' @param points A numeric vector representing points (indices of the vector x) to be labelled on the coverage plot.
#' @param distributions A PeakToolResult object. If no `points` are given, `points` will be extracted from `distributions`.
#' @param ... Additional parameters for the BoutrosLab.plotting.general function create.scatteplot
#' TODO: Add a parameter that allows for plotting distributions from something other than majority vote
#'
#'
#' @return Coverage plot, a Trellis object. For further details, see the 'Lattice' R package.
#' @export
#'
#' @examples
create.coverageplot = function(
  x,
  points = NULL,
  distributions = NULL,
  ...
){

  # Error checking
  # if(!is.null(distributions)), check that this is the right class

  # x
  if(!is.null(names(x))){
    labels.x = names(x)
    labels.x <- tryCatch(
        { labels.x = as.numeric(labels.x) },
        warning = function(cond) {
          message("Warning message:")
          message("Vector names are not coercible to numeric.")
          message("Using default indices.")
          return(1:length(x))
        }
      )
  } else {
    labels.x = 1:length(x)
  }
  plotting.data = data.frame("x" = x, "labels.x" = labels.x, "dist" = "coverage")

  # points
  if(!is.null(points)){
    points.x = labels.x[points]
    points.y = x[points]
  } else if( !is.null(distributions)){
    points = distributions[['p']]
    points = c(points[,'start'], points[,'end'])
    points.x = labels.x[points]
    points.y = x[points]
  } else {
    points.x = NULL
    points.y = NULL
  }

  # distribution
  # Extract fitted distributions of majority vote model
  if(!is.null(distributions)){
    models = distributions[['models']]
    majority_vote = lapply(models, `[[`, "majority.vote")
    distribution_plotting_data = lapply(majority_vote, function(m) {
      x = seq(m$seg.start, m$seg.end, by = 1)
      dens = m$dens(x = seq_along(x), mpar = m$par)
      data.frame("x" = dens, "labels.x" = labels.x[x], "dist" = m$dist)
    })
    distribution_plotting_data = do.call('rbind.data.frame', distribution_plotting_data)
    plotting.data = rbind(plotting.data, distribution_plotting_data)
  }

  # Factoring plotting data distribution
  plotting.data$dist = factor(plotting.data$dist, levels = c("coverage", "norm", "gamma", "unif"))

  # Plotting
  plt = BoutrosLab.plotting.general::create.scatterplot(
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
    ...
  )

  # Returning plotting object
  return(plt)
}
