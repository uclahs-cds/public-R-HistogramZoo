max.uniform.col <- 'firebrick3';
remove.low.entropy.col <- 'dodgerblue2';

common.legend <- legend.grob(
    list(
        legend = list(
            colours = c(max.uniform.col, remove.low.entropy.col),
            title = "Parameters",
            labels = c('Max uniform', 'Remove low entropy'),
            size = 3,
            title.cex = 1,
            label.cex = 1,
            border = 'black'
            ),
        legend = list(
            colours = distribution_colours[c('norm', 'gamma', 'unif')],
            title = "Distribution",
            labels = c("Normal","Gamma","Uniform"),
            size = 3,
            title.cex = 1,
            label.cex = 1,
            border = 'black'
            )
        ),
    label.cex = 1.5,
    title.cex = 1.5,
    size = 3,
    layout = c(1,2)
    )

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
    cov.data$max_uniform <- ifelse(cov.data$max_uniform, 'firebrick3', 'white');
  }
  if ('remove_low_entropy' %in% colnames(cov.data)) {
    cov.data$remove_low_entropy <- ifelse(cov.data$remove_low_entropy, 'dodgerblue', 'white');
  }
  if ('actual_dist' %in% colnames(cov.data)) {
    cov.data$actual_dist <- unname(distribution_colours[cov.data$actual_dist]);
  }
  if ('N_decile' %in% colnames(cov.data)) {
    cov.data$N_decile <- assign.cov.factor.col(cov.data$N_decile, c('white', 'slateblue4'));
  }
  if ('noise_decile' %in% colnames(cov.data)) {
    cov.data$noise_decile <- assign.cov.factor.col(cov.data$noise_decile, c('white', 'darkorange1'));
  }
  if ('eps_decile' %in% colnames(cov.data)) {
    cov.data$eps_decile <- assign.cov.factor.col(cov.data$eps_decile, c('white', 'seagreen3'));
  }
  if ('param_decile' %in% colnames(cov.data)) {
    cov.data$param_decile <- assign.cov.factor.col(cov.data$param_decile, c('white', 'yellowgreen'));
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
sim.plot.overall.accuracy <- function(x, ...) {
  # Accuracy by distribution
  overall.accuracy <- x[
      , .(
        accuracy_dist = mean(actual_dist == dist, na.rm = TRUE),
        accuracy_peaks = mean(num_segments == 1, na.rm = TRUE),
        N = mean(N),
        eps = mean(eps),
        noise = mean(noise)
        ), by = .(metric, actual_dist, max_uniform, remove_low_entropy)
      ];

  overall.accuracy.dist.wide <- dcast(
    overall.accuracy,
    max_uniform + remove_low_entropy + actual_dist ~ metric,
    value.var = 'accuracy_dist'
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

  create.multipanelplot(
    list(overall.accuracy.cov.heatmap, overall.accuracy.dist.heatmap),
    plot.objects.widths = c(0.1, 1),
    x.spacing = c(-0.25, 0),
    width = 12,
    legend = list(right = list(
      fun = common.legend
        )
      ),
    layout.width = 2,
    layout.height = 1,
    xlab.label = 'Distribution accuracy',
    ...
    )
}

sim.plot.quantile.accuracy <- function(x, ...) {
  decile.accuracy <- x[
    , .(
      accuracy_dist = mean(actual_dist == dist, na.rm = TRUE),
      accuracy_peaks = mean(num_segments == 1, na.rm = TRUE),
      N = .N
      ), by = .(metric, actual_dist, max_uniform, remove_low_entropy, N_decile, noise_decile, eps_decile)
    ]

  decile.accuracy <- decile.accuracy[!max_uniform & !remove_low_entropy, ];

  decile.wide.accuracy.dist <- dcast(
    decile.accuracy,
    actual_dist + max_uniform + remove_low_entropy + N_decile + noise_decile + eps_decile ~ metric,
    value.var = 'accuracy_dist'
    )

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

  create.multipanelplot(
    list(decile.accuracy.cov.heatmap, decile.wide.accuracy.dist.heatmap),
    plot.objects.widths = c(0.1, 1),
    x.spacing = c(-0.25, 0),
    width = 12,
    legend = list(right = list(
      fun = common.legend
        )
      ),
    layout.width = 2,
    layout.height = 1,
    xlab.label = 'Distribution accuracy',
    ...
    )
}
