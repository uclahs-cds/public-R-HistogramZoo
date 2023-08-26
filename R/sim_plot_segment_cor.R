#' Creates a Spearman correlation heatmap for metric and variables of interest including N, epslion, noise and parameter
#' Automatically separates by both distribution and `remove_low_entropy`
#'
#' @param x a data frame of simulation results with columns named for `cor.var` and simulation parameters including `N`, `eps`, `noise` and `param`
#' @param cluster whether or not to cluster the heatmap (using diana)
#' @param cor.var target variable for correlation, out of `jaccard` and `prob_segment`
#' @param ... additional arguments to `BoutrosLab.plotting.general::create.multipanelplot`
#'
#' @export
sim.plot.segment.cor <- function(
    x,
    cluster = FALSE,
    cor.var = c('jaccard', 'prob_segment'),
    ...
){
  cor.var <- match.arg(cor.var)

  # Calculating spearman's correlation
  seg.cor.data <- as.data.frame(x[
    ,
    .(
      N_cor = cor(N, get(cor.var), method = 'spearman'),
      eps_cor = cor(eps, get(cor.var), method = 'spearman'),
      noise_cor = cor(noise, get(cor.var), method = 'spearman'),
      param_cor = cor(param, get(cor.var), method = 'spearman')
      )
    , by = .(remove_low_entropy, actual_dist)
    ][
      order(actual_dist, remove_low_entropy)
    ])

  cor.cols <- colnames(seg.cor.data)[grepl('_cor', colnames(seg.cor.data))]

  # Clustering and extracting variable order
  if (cluster) {
    seg.cor.order <- diana(seg.cor.data[, cor.cols])$order
  } else {
    seg.cor.order <- seq_len(nrow(seg.cor.data))
  }

  # Creating covariate heatmap
  cor.cov.heatmap <- sim.plot.heatmap.cov(
      seg.cor.data[seg.cor.order, c(
        'actual_dist', 'remove_low_entropy'
        )]
      );

  # Reordering dataframe
  cor.data <- seg.cor.data[seg.cor.order, cor.cols]

  # Adding text to the heatmap
  row.col.text <- cbind.data.frame(
    which(abs(cor.data) > 0.1, arr.ind = TRUE),
    value = round(cor.data[abs(cor.data) > 0.1], digits = 2)
    )
  row.col.text$row <- nrow(cor.data) - row.col.text$row + 1

  # X axis labels
  xaxis.lab <- c(
    expression('\u03c1'['N']),
    expression('\u03c1'['\u03B5']),
    expression('\u03c1'['noise']),
    expression('\u03c1'['param'])
    )

  # Correlation heatmap
  cov.hm <- create.heatmap(
    x = cor.data,
    same.as.matrix = TRUE,
    colourkey.cex = 2,
    clustering.method = 'none',
    row.pos = row.col.text$row,
    col.pos = row.col.text$col,
    cell.text = row.col.text$value,
    text.use.grid.coordinates = TRUE,
    xaxis.lab = xaxis.lab,
    xaxis.rot = 0,
    at = seq(-1, 1, length.out = 10),
    xaxis.tck = 0,
    yaxis.tck = 0,
    fill.colour = 'lightgrey'
    )

  # Multipanelplot
  create.multipanelplot(
    list(cor.cov.heatmap, cov.hm),
    plot.objects.widths = c(0.1, 1),
    x.spacing = c(-0.25, 0),
    main.cex = 2,
    layout.width = 2,
    layout.height = 1,
    xlab.label = 'Spearman correlation',
    legend = list(
        right = list(
          fun = common.sim.legend(
            include.legends = c('params', 'distributions'),
            params.to.include = c('remove_low_entropy'),
            cont.params.to.include = NULL
            )
          )
        ),
    ...
    )
}
