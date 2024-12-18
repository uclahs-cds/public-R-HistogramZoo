#' Generates a confusion matrix for predicted and true distribution
#'
#' @param pred a vector of predicted distributions
#' @param actual a vector of true distributions
#' @param eval_metrics evaluation metrics for confusion matrix generation
#'
#' @return a data table with columns `Eval Metric`, `dist` and `value` describing the evaluation metric for each distribution
confusion_matrix_dist <- function(
    pred,
    actual,
    eval_metrics = c(
      'Sensitivity', 'Specificity', 'Pos Pred Value', 'Neg Pred Value',
      'Precision', 'Recall', 'F1', 'Prevalence', 'Detection Rate',
      'Detection Prevalence', 'Balanced Accuracy'
      )
){
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

#' Creates an accuracy heatmap relative to metrics
#'
#' @param x a data frame of simulation results
#' @param acc accuracy variable; one of `dist` (distribution), `peaks` (number of peaks), `both`, `count`, and options from `caret::confusionMatrix`
#' @param cluster whether or not to cluster (using diana)
#' @param print.colour.key whether or not to print the colour key for the heatmap
#' @param group_vars variables to group in the covariate
#' @param legend a `legend.grob` object from `BoutrosLab.plotting.general`
#' @param xlab.label x axis label
#' @param ... additional parameters to be passed to `BoutrosLab.plotting.general::create.multipanelplot`
#'
#' @import data.table
#' @export
sim.plot.distribution.accuracy <- function(
    x,
    acc = c(
      'dist', 'peaks', 'both', 'count',
      # Confusion matrix metrics
      'Sensitivity', 'Specificity', 'Pos Pred Value', 'Neg Pred Value',
      'Precision', 'Recall', 'F1', 'Prevalence', 'Detection Rate',
      'Detection Prevalence', 'Balanced Accuracy'
      ),
    cluster = FALSE,
    print.colour.key = TRUE,
    group_vars = c(
      'actual_dist', 'max_uniform', 'remove_low_entropy',
      'jaccard_decile', 'N_decile', 'noise_decile', 'eps_decile',
      'interference_decile', 'proportion_decile'
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
    ...
){
  acc <- match.arg(acc);
  group_vars <- match.arg(group_vars, several.ok = TRUE); # Restricting the set of group vars for now
  # TODO: Error checking:
  # 1. Make sure all required columns are in data frame
  # 2. "actual_dist" *has to* be a covariate for caret functions
  # 3. NA's in the column of interest will lead to nonsensical results

  # Confusion matrix accuracy calculation
  is.cm.acc <- acc %in% c('Sensitivity', 'Specificity', 'Pos Pred Value', 'Neg Pred Value',
      'Precision', 'Recall', 'F1', 'Prevalence', 'Detection Rate',
      'Detection Prevalence', 'Balanced Accuracy')
  if (is.cm.acc) {
    confusion_matrix_vars <- setdiff(group_vars, c('actual_dist', 'dist'))
    decile.accuracy <- x[
      ,
      confusion_matrix_dist(
          dist,
          actual_dist,
          eval_metrics = ..acc
          ),
      by = c('metric', confusion_matrix_vars)
      ]

    group_vars[group_vars == "actual_dist"] <- "dist" # Accounting for confusion_matrix_dist output format
    value.var <- 'value'
  } else {

    # Accuracy of distribution and number of peaks
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

  # Generating a long-wide format of the data
  metrics <- unique(decile.accuracy$metric)
  metrics <- metrics[match(metrics, metric_ref[metric_ref %in% metrics])]

  long.wide.formula <- paste(paste0(group_vars, collapse = ' + '), 'metric', sep = ' ~ ')
  decile.wide.accuracy.dist <- dcast(
    data = decile.accuracy,
    formula = long.wide.formula,
    value.var = value.var
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

  # accuracy heatmap
  decile.wide.accuracy.dist.heatmap <- create.heatmap(
    decile.wide.accuracy.dist[diana.acc.clust, metrics],
    same.as.matrix = TRUE,
    clustering.method = 'none',
    colour.scheme = if (acc != 'count') acc.colour.scheme else c('white', 'forestgreen'),
    xaxis.lab = metrics,
    xaxis.rot = 90,
    at = if (acc == 'count') NULL else seq(0, 1, length.out = 20),
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
    list(decile.accuracy.cov.heatmap, decile.wide.accuracy.dist.heatmap),
    legend = legend,
    main.cex = 2,
    layout.width = 2,
    layout.height = 1,
    xlab.label = xlab.label,
    ...
    )
}
