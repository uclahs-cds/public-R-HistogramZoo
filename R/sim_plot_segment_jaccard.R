jaccard.colour.scheme <- c('white', 'dodgerblue2');

#' @export
sim.plot.segment.eval <- function(
    x,
    cluster = FALSE,
    print.colour.key = TRUE,
    group_vars = c(
      'remove_low_entropy',
      'N_decile', 'noise_decile',
      'eps_decile'
      ),
    target = c('median_jaccard', 'count', 'prob_segment'),
    legend = list(
      right = list(
        fun = common.sim.legend(
          include.legends = c('params', 'distributions', 'quantiles'),
          params.to.include = intersect(group_vars, c('remove_low_entropy')),
          cont.params.to.include = group_vars[grepl('_decile', group_vars)]
          )
        )
      ),
    ...
    ) {
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


  jaccard.cov.heatmap <- sim.plot.heatmap.cov(
      data.frame(res.wide[diana.jaccard.clust, ..group_vars])
      );

  heatmap.at <- if (target == 'count') NULL else seq(0, 1, length.out = 20);

  jaccard.heatmap <- create.heatmap(
      res.wide[diana.jaccard.clust, c('gamma', 'norm', 'unif')],
      same.as.matrix = TRUE,
      print.colour.key = print.colour.key,
      clustering.method = 'none',
      colour.scheme = jaccard.colour.scheme,
      xaxis.lab = c('gamma', 'norm', 'unif'),
      xaxis.rot = 45,
      at = heatmap.at,
      xaxis.tck = 0,
      yaxis.tck = 0,
      fill.colour = 'lightgrey'
      );

  create.multipanelplot(
      list(jaccard.cov.heatmap, jaccard.heatmap),
      plot.objects.widths = c(0.1, 1),
      x.spacing = c(-0.25, 0),
      width = 12,
      legend = legend,
      main.cex = 2,
      main = if (cluster) 'DIANA Clustered' else 'No Clustering',
      layout.width = 2,
      layout.height = 1,
      ...
      )
  }
