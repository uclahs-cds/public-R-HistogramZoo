covariate.col <- list(
  max.uniform = 'firebrick3',
  remove.low.entropy = 'dodgerblue2',
  N.decile = 'slateblue4',
  noise.decile = 'darkorange1',
  eps.decile = 'seagreen3',
  param.decile = 'yellowgreen'
  )


common.sim.legend <- function(
    include.legends = c('params', 'distributions', 'quantiles'),
    params.to.include = c(
      'max_uniform', 'remove_low_entropy'
      ),
    cont.params.to.include = c(
      'N_decile', 'noise_decile',
      'eps_decile', 'param_decile'
    )
  ) {
  include.legends <- match.arg(include.legends, several.ok = TRUE);
  params.to.include <- match.arg(params.to.include, several.ok = TRUE);
  if (! is.null(cont.params.to.include)) {
    cont.params.to.include <- match.arg(cont.params.to.include, several.ok = TRUE);
    }

  legend.list <- list()
  if ('params' %in% include.legends) {
    params.select <- match(params.to.include, c('max_uniform', 'remove_low_entropy'));

    legend.list <- c(
      legend.list,
      list(
        legend = list(
            colours = c(covariate.col$max.uniform, covariate.col$remove.low.entropy)[params.select],
            title = "Parameters",
            labels = c('Max uniform', 'Remove low entropy')[params.select],
            size = 3,
            title.cex = 1,
            label.cex = 1,
            border = 'black'
            )
        )
      )
    }
  if ('distributions' %in% include.legends) {
    legend.list <- c(
      legend.list,
      list(
        legend = list(
            colours = distribution_colours[c('norm', 'gamma', 'unif')],
            title = "Distribution",
            labels = c("Normal","Gamma","Uniform"),
            size = 3,
            title.cex = 1,
            label.cex = 1,
            border = 'black'
            )
        )
      )
    }

  if (! is.null(cont.params.to.include)) {
    cont.params <- gsub('_', '.', cont.params.to.include)

    cont.params.legend <- lapply(cont.params, function(p) {
      list(
        colours = c('white', covariate.col[[p]]),
        labels = c('0', '1'), # TODO
        title = sub('[.]decile', '', p),
        angle = -90,
        width = 6,
        continuous = TRUE
        )
      })

    names(cont.params.legend) <- rep('legend', length(cont.params.legend))

    legend.list <- c(legend.list, cont.params.legend)
    }

  BoutrosLab.plotting.general::legend.grob(
    legends = legend.list,
    label.cex = 1.5,
    title.cex = 1.5,
    size = 3,
    title.just = 'left',
    layout = c(1,length(legend.list))
    )
  }



assign.cov.factor.col <- function(x, ramp.colors) {
  x.levels <- levels(x);
  k <- length(x.levels);

  colormap <- setNames(
    colorRampPalette(ramp.colors)(k),
    x.levels
    );

  colormap[x];
}

sim.plot.heatmap.cov <- function(cov.data) {
  if ('max_uniform' %in% colnames(cov.data)) {
    cov.data$max_uniform <- ifelse(cov.data$max_uniform, covariate.col$max.uniform, 'white');
  }
  if ('remove_low_entropy' %in% colnames(cov.data)) {
    cov.data$remove_low_entropy <- ifelse(cov.data$remove_low_entropy, covariate.col$remove.low.entropy, 'white');
  }
  if ('actual_dist' %in% colnames(cov.data)) {
    cov.data$actual_dist <- unname(distribution_colours[cov.data$actual_dist]);
  }
  if ('N_decile' %in% colnames(cov.data)) {
    cov.data$N_decile <- assign.cov.factor.col(cov.data$N_decile, c('white', covariate.col$N.decile));
  }
  if ('noise_decile' %in% colnames(cov.data)) {
    cov.data$noise_decile <- assign.cov.factor.col(cov.data$noise_decile, c('white', covariate.col$noise.decile));
  }
  if ('eps_decile' %in% colnames(cov.data)) {
    cov.data$eps_decile <- assign.cov.factor.col(cov.data$eps_decile, c('white', covariate.col$eps.decile));
  }
  if ('param_decile' %in% colnames(cov.data)) {
    cov.data$param_decile <- assign.cov.factor.col(cov.data$param_decile, c('white', covariate.col$param.decile));
  }

  if (ncol(cov.data) == 1) cov.data <- cbind(cov.data[, 1], cov.data[, 1]);

  create.heatmap(
    cov.data,
    same.as.matrix = TRUE,
    input.colours = TRUE,
    clustering.method = 'none',
    print.colour.key = FALSE,
    xaxis.tck = 0,
    yaxis.tck = 0
    )
}

#' @export
#' @import data.table
sim.plot.overall.accuracy <- function(
    x,
    acc = c('dist', 'peaks', 'both'),
    ...) {
  acc <- match.arg(acc);
  # Accuracy by distribution
  overall.accuracy <- x[
      , .(
        accuracy_dist = mean(actual_dist == dist, na.rm = TRUE),
        accuracy_peaks = mean(num_segments == 1, na.rm = TRUE),
        accuracy_both = mean(actual_dist == dist & num_segments == 1, na.rm = TRUE),
        N = mean(N),
        eps = mean(eps),
        noise = mean(noise)
        ), by = .(metric, actual_dist, max_uniform, remove_low_entropy)
      ];

  overall.accuracy.dist.wide <- dcast(
    overall.accuracy,
    max_uniform + remove_low_entropy + actual_dist ~ metric,
    value.var = paste0('accuracy_', acc)
    )

  overall.accuracy.cov.heatmap <- sim.plot.heatmap.cov(
    data.frame(overall.accuracy.dist.wide[, c('max_uniform', 'remove_low_entropy', 'actual_dist')])
    );

  overall.accuracy.dist.heatmap <- create.heatmap(
    overall.accuracy.dist.wide[, ..metrics],
    same.as.matrix = TRUE,
    clustering.method = 'none',
    colour.scheme = c('white', 'red'),
    xaxis.lab = metrics,
    xaxis.rot = 45,
    at = seq(0, 1, length.out = 20),
    xaxis.tck = 0,
    yaxis.tck = 0,
    fill.colour = 'lightgrey'
    );

  xlab.label <- switch (acc,
    'dist' = 'Distribution accuracy',
    'peaks' = 'Number of peaks accuracy',
    'both' = 'Correct peak and distribution accuracy'
    )

  create.multipanelplot(
    list(overall.accuracy.cov.heatmap, overall.accuracy.dist.heatmap),
    plot.objects.widths = c(0.1, 1),
    x.spacing = c(-0.25, 0),
    width = 12,
    legend = list(right = list(
      fun = common.sim.legend(
          cont.params.to.include = NULL
          )
        )
      ),
    layout.width = 2,
    layout.height = 1,
    xlab.label = xlab.label,
    ...
    )
}

sim.plot.quantile.accuracy <- function(
    x,
    acc = c('dist', 'peaks', 'both'),
    sort.cols = NULL, ...) {
  acc <- match.arg(acc);
  decile.accuracy <- x[
    , .(
      accuracy_dist = mean(actual_dist == dist, na.rm = TRUE),
      accuracy_peaks = mean(num_segments == 1, na.rm = TRUE),
      accuracy_both = mean(actual_dist == dist & num_segments == 1, na.rm = TRUE),
      N = .N
      ), by = .(metric, actual_dist, max_uniform, remove_low_entropy, N_decile, noise_decile, eps_decile)
    ]

  decile.accuracy <- decile.accuracy[!max_uniform & !remove_low_entropy, ];

  decile.wide.accuracy.dist <- dcast(
    decile.accuracy,
    actual_dist + max_uniform + remove_low_entropy + N_decile + noise_decile + eps_decile ~ metric,
    value.var = paste0('accuracy_', acc)
    )

  cont.cols <- colnames(decile.wide.accuracy.dist)[grepl('_decile', colnames(decile.wide.accuracy.dist))];

  if (!is.null(sort.cols)) {
    setorderv(decile.wide.accuracy.dist, sort.cols);
  }

  decile.accuracy.cov.heatmap <- sim.plot.heatmap.cov(
    #data.frame(decile.wide.accuracy.dist[, c('max_uniform', 'remove_low_entropy', 'actual_dist'), drop = FALSE])
    data.frame(decile.wide.accuracy.dist[, c('actual_dist', 'N_decile', 'noise_decile', 'eps_decile')])
    );

  decile.wide.accuracy.dist.heatmap <- create.heatmap(
    decile.wide.accuracy.dist[, ..metrics],
    same.as.matrix = TRUE,
    clustering.method = 'none',
    colour.scheme = c('white', 'red'),
    xaxis.lab = metrics,
    xaxis.rot = 45,
    at = seq(0, 1, length.out = 20),
    xaxis.tck = 0,
    yaxis.tck = 0,
    fill.colour = 'lightgrey'
    );

  xlab.label <- switch (acc,
    'dist' = 'Distribution accuracy',
    'peaks' = 'Number of peaks accuracy',
    'both' = 'Correct peak and distribution accuracy'
    )

  create.multipanelplot(
    list(decile.accuracy.cov.heatmap, decile.wide.accuracy.dist.heatmap),
    plot.objects.widths = c(0.1, 1),
    x.spacing = c(-0.25, 0),
    width = 12,
    legend = list(right = list(
      fun = common.sim.legend(
        include.legends = c('distributions', 'quantiles'),
        cont.params.to.include = cont.cols
        )
      )
    ),
    layout.width = 2,
    layout.height = 1,
    xlab.label = xlab.label,
    ...
    )
}
