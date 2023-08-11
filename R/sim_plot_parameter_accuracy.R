

#' Creates a parameter accuracy heatmap relative to metrics
#'
#' @param x a data frame of simulation results
#' @param acc regression accuracy metric from the package `Metrics`; one of `mse` (mean squared error), `mdae` (median absolute error), `rmse` (root mean squared error) and `mae` (mean absolute error)
#' @param cluster whether or not to cluster (using diana)
#' @param print.colour.key whether or not to print the colour key for the heatmap
#' @param group_vars variables to group in the covariate
#' @param legend a `legend.grob` object from `BoutrosLab.plotting.general`
#' @param xlab.label x axis label
#' @param ... additional parameters to be passed to `BoutrosLab.plotting.general::create.multipanelplot`
#'
#' @import data.table
#' @import Metrics
#' @export
sim.plot.parameter.accuracy <- function(
    x,
    acc = c('mse', 'mdae', 'rmse', 'mae'),
    cluster = FALSE,
    print.colour.key = TRUE,
    group_vars = c(
      'actual_dist', 'max_uniform', 'remove_low_entropy',
      'jaccard_decile', 'N_decile', 'noise_decile', 'eps_decile'
    ),
    legend = list(
      right = list(
        fun = common.sim.legend(
          include.legends = c('params', 'distributions', 'quantiles'),
          params.to.include = intersect(group_vars, c('max_uniform', 'remove_low_entropy')),
          cont.params.to.include = group_vars[grepl('_decile', group_vars)]
        )
      )
    ),
    xlab.label = switch(
      acc,
      "mse" = "Mean Squared Error",
      "mae" = "Mean Average Error",
      "mdae" = "Median Absolute Error",
      "rmse" = "Root Mean Squared Error"
    ),
    ...
){
  acc <- match.arg(acc);
  group_vars <- match.arg(group_vars, several.ok = TRUE); # Restricting the set of group vars for now
  # TODO: Error checking:
  # 1. Make sure all required columns are in data frame
  # 2. "actual_dist" *has to* be a covariate for caret functions
  
  # Accuracy of parameter estimate
  decile.accuracy <- x[
    , 
    .(
      "value" = do.call(acc, list(param, param_fit))
     ), by = c('metric', group_vars)
  ]

  # Generating a long-wide format of the data
  metrics <- unique(decile.accuracy$metric)
  
  long.wide.formula <- paste(paste0(group_vars, collapse = ' + '), 'metric', sep = ' ~ ')
  decile.wide.accuracy.dist <- dcast(
    data = decile.accuracy,
    formula = long.wide.formula,
    value.var = "value"
  )
  decile.wide.accuracy.dist <- as.data.frame(decile.wide.accuracy.dist)
  
  # Clustering data
  if (cluster) {
    cluster.mat <- decile.wide.accuracy.dist[, metrics]
    na.cells <- is.na(decile.wide.accuracy.dist[, metrics])
    cluster.mat[na.cells] <- -1
    
    diana.acc.clust <- diana(
      cluster.mat
    )$order
  } else {
    diana.acc.clust <- seq_len(nrow(decile.wide.accuracy.dist))
  }
  
  # covariate heatmap
  decile.accuracy.cov.heatmap <- sim.plot.heatmap.cov(
    decile.wide.accuracy.dist[diana.acc.clust, group_vars]
  )
  
  colour_bar <- range(decile.accuracy$value)
  colour_bar <- c(floor(colour_bar[1]), ceiling(colour_bar[2]))
  
  # accuracy heatmap
  decile.wide.accuracy.dist.heatmap <- create.heatmap(
    decile.wide.accuracy.dist[diana.acc.clust, metrics],
    same.as.matrix = TRUE,
    clustering.method = 'none',
    print.colour.key = print.colour.key,
    colour.scheme = c('white', 'darkmagenta'),
    xaxis.lab = metrics,
    xaxis.rot = 90,
    at = seq(colour_bar[1], colour_bar[2], length.out = 20),
    xaxis.tck = 0,
    yaxis.tck = 0,
    fill.colour = 'lightgrey'
  );
  
  create.multipanelplot(
    list(decile.accuracy.cov.heatmap, decile.wide.accuracy.dist.heatmap),
    plot.objects.widths = c(0.3, 1),
    x.spacing = c(-0.25, 0),
    legend = legend,
    main.cex = 2,
    layout.width = 2,
    layout.height = 1,
    xlab.label = xlab.label,
    ...
  )
}
