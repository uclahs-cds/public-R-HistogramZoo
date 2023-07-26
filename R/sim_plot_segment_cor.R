#' @export
sim.plot.segment.cor <- function(
    x,
    cluster = FALSE,
    cor.var = c('jaccard', 'prob_segment'),
    ...
    ) {
  cor.var <- match.arg(cor.var)

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

  if (cluster) {
    seg.cor.order <- diana(seg.cor.data[, cor.cols])$order
  } else {
    seg.cor.order <- seq_len(nrow(seg.cor.data))
  }


  cor.cov.heatmap <- sim.plot.heatmap.cov(
      seg.cor.data[seg.cor.order, c(
        'actual_dist', 'remove_low_entropy'
        )]
      );

  cor.data <- seg.cor.data[seg.cor.order, cor.cols]

  row.col.text <- cbind.data.frame(
    which(abs(cor.data) > 0.1, arr.ind = TRUE),
    value = round(cor.data[abs(cor.data) > 0.1], digits = 2)
    )
  row.col.text$row <- nrow(cor.data) - row.col.text$row + 1

  xaxis.lab <- c(
    expression('\u03c1'['N']),
    expression('\u03c1'['\u03B5']),
    expression('\u03c1'['noise']),
    expression('\u03c1'['param'])
    )

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

  create.multipanelplot(
    list(cor.cov.heatmap, cov.hm),
    plot.objects.widths = c(0.1, 1),
    x.spacing = c(-0.25, 0),
    width = 12,
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
