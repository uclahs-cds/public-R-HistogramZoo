jaccard.colour.scheme <- c('white', 'royalblue4');


#' Creates a heatmap of a metric vs. distributions for simulation data split by different parameters
#'
#' @param x simulation dataset in the form of a data frame with corresponding column names for `actual_dist` (distributions), `group_vars`and `target`
#' @param cluster whether or not to cluster using the method `diana`
#' @param print.colour.key whether or not to print the colour key
#' @param group_vars columns used to group variable of interest
#' @param target column name of variable of interest to aggregate and visualize
#' @param legend a `legend.grob` from `BoutrosLab.plotting.general`
#' @param ... additional parameters to be passed to `create.multipanelplot`
#'
#' @export
sim.plot.segment.eval <- function(
    x,
    cluster = FALSE,
    print.colour.key = TRUE,
    group_vars = c(
      'remove_low_entropy', 'max_uniform',
      'N_decile', 'noise_decile',
      'eps_decile'
      ),
    target = c('median_jaccard', 'count', 'prob_segment'),
    legend = list(
      right = list(
        fun = common.sim.legend(
          include.legends = c('params', 'distributions', 'quantiles'),
          params.to.include = intersect(group_vars, c('remove_low_entropy', 'max_uniform')),
          cont.params.to.include = group_vars[grepl('_decile', group_vars)],
          simulation.params = list(
            'N' = c(25, 500),
            'eps' = c(0.5, 2),
            'noise' = c(.05, .5)
            )
          )
        )
      ),
    ...
){

  target <- match.arg(target);
  if (target == 'median_jaccard') {
      res <- x[
      ,
      .(
        median_metric = median(jaccard),
        count = .N
        ), by = c('actual_dist', group_vars)
      ]
      target <- 'median_metric'
  } else if (target == 'prob_segment') {
      res <- x[
      ,
      .(
        median_metric = median(prob_segment),
        count = .N
        ), by = c('actual_dist', group_vars)
      ]
      target <- 'median_metric'
  } else {
      res <- x[
      ,
      .(
        count = .N
        ), by = c('actual_dist', group_vars)
      ]
    }


  long.wide.formula <- paste(paste0(group_vars, collapse = ' + '), 'actual_dist', sep = ' ~ ')
  res.wide <- dcast(
      data = res,
      formula = long.wide.formula,
      value.var = target
      )

  # Scale
  if (target == 'count') {
    res.wide[, c('gamma', 'norm', 'unif')] <- as.data.frame(scale(res.wide[, c('gamma', 'norm', 'unif')]))
    }

  if (cluster) {
      diana.jaccard.clust <- diana(
        res.wide[, c('gamma', 'norm', 'unif')]
      )$order
  } else {
      diana.jaccard.clust <- seq_len(nrow(res.wide))
    }

  # Generating covariate heatmap
  jaccard.cov.heatmap <- sim.plot.heatmap.cov(
      data.frame(res.wide[diana.jaccard.clust, ..group_vars])
      );

  heatmap.at <- if (target == 'count') NULL else seq(0, 1, length.out = 20);

  jaccard.heatmap <- create.heatmap(
      res.wide[diana.jaccard.clust, c('gamma', 'norm', 'unif')],
      same.as.matrix = TRUE,
      clustering.method = 'none',
      colour.scheme = jaccard.colour.scheme,
      xaxis.lab = distribution_names[c('gamma', 'norm', 'unif')],
      xaxis.rot = 90,
      at = heatmap.at,
      xaxis.tck = 0,
      yaxis.tck = 0,
      fill.colour = 'lightgrey',
      # Colourkey
      print.colour.key = print.colour.key,
      # colourkey.labels.at = NULL,
      # colourkey.labels = NULL,
      colourkey.cex = 1.5
      );

  create.multipanelplot(
      list(jaccard.cov.heatmap, jaccard.heatmap),
      legend = legend,
      main.cex = 2,
      main = if (cluster) 'DIANA Clustered' else 'No Clustering',
      layout.width = 2,
      layout.height = 1,
      ...
      )
  }
