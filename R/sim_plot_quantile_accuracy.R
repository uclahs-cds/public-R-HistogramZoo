confusion_matrix_dist <- function(
    pred,
    actual,
    eval_metrics = c(
      'Sensitivity', 'Specificity', 'Pos Pred Value', 'Neg Pred Value',
      'Precision', 'Recall', 'F1', 'Prevalence', 'Detection Rate',
      'Detection Prevalence', 'Balanced Accuracy'
      )
    ) {
  eval_metrics <- match.arg(eval_metrics, several.ok = TRUE)
#  if (length(actual) == 1) actual <- rep(actual, length(pred))
  res <- as.data.table(
    t(caret::confusionMatrix(
      factor(pred, levels = c('gamma', 'norm', 'unif', 'gamma_flip')),
      factor(actual, levels = c('gamma', 'norm', 'unif', 'gamma_flip'))
      )$byClass[-4, eval_metrics]),
    keep.rownames = TRUE
    )
  colnames(res) <- sub('Class: ', '', colnames(res))
  names(res)[names(res) == 'rn'] <- 'Eval Metric'
  melt(
    res,
    id.vars = setdiff(colnames(res), c('gamma', 'norm', 'unif')),
    variable.name = 'dist'
    )
}

sim.plot.quantile.cm.metric <- function(
    x,
    eval_metric = c(
      'Sensitivity', 'Specificity', 'Pos Pred Value', 'Neg Pred Value',
      'Precision', 'Recall', 'F1', 'Prevalence', 'Detection Rate',
      'Detection Prevalence', 'Balanced Accuracy'
      ),
    cluster = FALSE,
    print.colour.key = TRUE,
    group_vars = c(
      'max_uniform', 'remove_low_entropy',
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
    ...
  ) {
  eval_metric <- match.arg(eval_metric)
  group_vars <- setdiff(group_vars, c('actual_dist', 'dist'))
  cm_matrix <- x[
    ,
    confusion_matrix_dist(
        dist,
        actual_dist,
        eval_metrics = ..eval_metric
        ),
    by = c('metric', group_vars)
    ]

  long.wide.formula <- paste(paste0(c(group_vars, 'dist'), collapse = ' + '), 'metric', sep = ' ~ ')
  cm_matrix_wide <- dcast(
    cm_matrix,
    long.wide.formula,
    value.var = 'value'
    )


  browser()
  }

#' @import data.table
#' @export
sim.plot.quantile.accuracy <- function(
    x,
    acc = c(
      'dist', 'peaks', 'both', 'count',
      # Confusion matrix metrics
      'Sensitivity', 'Specificity', 'Pos Pred Value', 'Neg Pred Value',
      'Precision', 'Recall', 'F1', 'Prevalence', 'Detection Rate',
      'Detection Prevalence', 'Balanced Accuracy'
      ),
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
      'both' = 'Correct peak and distribution accuracy',
        {
        acc
        }
      ),
    ...) {
  acc <- match.arg(acc);

  is.cm.acc <- acc %in% c('Sensitivity', 'Specificity', 'Pos Pred Value', 'Neg Pred Value',
      'Precision', 'Recall', 'F1', 'Prevalence', 'Detection Rate',
      'Detection Prevalence', 'Balanced Accuracy')
  if (is.cm.acc) {
    group_vars <- setdiff(group_vars, c('actual_dist', 'dist'))
    decile.accuracy <- x[
      ,
      confusion_matrix_dist(
          dist,
          actual_dist,
          eval_metrics = ..acc
          ),
      by = c('metric', group_vars)
      ]

    # Correct ordering of covariates
    group_vars <- c(
      group_vars[group_vars %in% c('max_uniform', 'remove_low_entropy')],
      'dist',
      group_vars[! group_vars %in% c('max_uniform', 'remove_low_entropy')]
      )
    value.var <- 'value'
  } else {
    decile.accuracy <- x[
    , .(
      accuracy_dist = mean(actual_dist == dist, na.rm = TRUE),
      accuracy_peaks = mean(num_segments == 1, na.rm = TRUE),
      accuracy_both = mean(actual_dist == dist & num_segments == 1, na.rm = TRUE),
      count = .N
      ), by = c('metric', group_vars)
    ]

    value.var <- if (acc == 'count') acc else paste0('accuracy_', acc)
  }

  metrics <- unique(decile.accuracy$metric)

  long.wide.formula <- paste(paste0(group_vars, collapse = ' + '), 'metric', sep = ' ~ ')
  decile.wide.accuracy.dist <- dcast(
    decile.accuracy,
    long.wide.formula,
    value.var = value.var
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
    colour.scheme = if (acc != 'count') acc.colour.scheme else c('white', 'forestgreen'),
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
