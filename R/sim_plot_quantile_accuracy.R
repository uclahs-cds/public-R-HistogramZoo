#' @import data.table
#' @export
sim.plot.quantile.accuracy <- function(
    x,
    acc = c('dist', 'peaks', 'both', 'count'),
    sort.cols = NULL,
    cluster = FALSE,
    print.colour.key = TRUE,
    group_vars = c(
      'actual_dist', 'max_uniform', 'remove_low_entropy',
      'jaccard_decile', 'N_decile', 'noise_decile'
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
    xlab.label = switch (acc,
      'dist' = 'Distribution accuracy',
      'peaks' = 'Number of peaks accuracy',
      'both' = 'Correct peak and distribution accuracy'
      ),
    ...) {
  acc <- match.arg(acc);
  decile.accuracy <- x[
    , .(
      accuracy_dist = mean(actual_dist == dist, na.rm = TRUE),
      accuracy_peaks = mean(num_segments == 1, na.rm = TRUE),
      accuracy_both = mean(actual_dist == dist & num_segments == 1, na.rm = TRUE),
      count = .N
      ), by = c('metric', group_vars)
    ]

  metrics <- unique(decile.accuracy$metric)

  # decile.accuracy <- decile.accuracy[(!max_uniform), ];

  long.wide.formula <- paste(paste0(group_vars, collapse = ' + '), 'metric', sep = ' ~ ')
  decile.wide.accuracy.dist <- dcast(
    decile.accuracy,
    long.wide.formula,
    value.var = if (acc == 'count') acc else paste0('accuracy_', acc)
    )
  decile.wide.accuracy.dist <- as.data.frame(decile.wide.accuracy.dist)

  cont.cols <- colnames(decile.wide.accuracy.dist)[grepl('_decile', colnames(decile.wide.accuracy.dist))];
  if (cluster) {
    cluster.mat <- decile.wide.accuracy.dist[, metrics]
    na.cells <- is.na(decile.wide.accuracy.dist[, metrics])
    cluster.mat[na.cells] <- -1

    # diana.acc.clust <- hclust(
    #   dist(
    #     cluster.mat
    #   )
    # )$order
    diana.acc.clust <- diana(
      cluster.mat
      )$order
  } else {
      diana.acc.clust <- seq_len(nrow(decile.wide.accuracy.dist))
    }

  decile.accuracy.cov.heatmap <- sim.plot.heatmap.cov(
    decile.wide.accuracy.dist[diana.acc.clust, group_vars]
    )

  decile.wide.accuracy.dist.heatmap <- create.heatmap(
    decile.wide.accuracy.dist[diana.acc.clust, metrics],
    same.as.matrix = TRUE,
    clustering.method = 'none',
    print.colour.key = print.colour.key,
    colour.scheme = acc.colour.scheme,
    xaxis.lab = metrics,
    xaxis.rot = 45,
    at = if (acc == 'count') NULL else seq(0, 1, length.out = 20),
    xaxis.tck = 0,
    yaxis.tck = 0,
    fill.colour = 'lightgrey'
    );

  create.multipanelplot(
    list(decile.accuracy.cov.heatmap, decile.wide.accuracy.dist.heatmap),
    plot.objects.widths = c(0.1, 1),
    x.spacing = c(-0.25, 0),
    width = 12,
    legend = legend,
    main.cex = 2,
    layout.width = 2,
    layout.height = 1,
    xlab.label = xlab.label,
    ...
    )
}
