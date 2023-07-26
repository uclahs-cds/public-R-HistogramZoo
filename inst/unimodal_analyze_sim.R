library(HistogramZoo);
library(data.table);
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

base.path <- HistogramZoo:::load.config()$root.path;
results.folder <- file.path(base.path, 'results');
merged.folder <- file.path(results.folder, 'merged_sims');
plots.folder <- file.path(base.path, 'plots');

sim.version <- c('v5', 'v6');
sim.merge.date <- '2023-07-24'

sim.version.regex <- paste0('[', paste0(sim.version, collapse = '|'), ']')
sim.version <- paste0(sim.version, collapse = '-')
metric.files <- list.files(
  path = merged.folder,
  pattern = paste0(sim.merge.date, '_.*', sim.version.regex, '.*noise.tsv'),
  full.names = TRUE
  )

unimodal.sim.metrics <- rbindlist(
  lapply(metric.files, data.table::fread, sep = '\t'),
  fill = TRUE
  )

mle.files <- list.files(
  path = merged.folder,
  pattern = paste0(sim.merge.date, '_.*', sim.version.regex, '.*mle'),
  full.names = TRUE
  )

unimodal.sim.metrics.mle <- rbindlist(
  lapply(mle.files, data.table::fread, sep = '\t'),
  fill = TRUE
  )

unimodal.sim.metrics.mle$max_uniform <- FALSE;

unimodal.sim <- rbindlist(
  list(
    unimodal.sim.metrics,
    unimodal.sim.metrics.mle
    ),
  fill = TRUE
  )

metrics <- c('mle', 'chisq', 'intersection', 'jaccard', 'ks', 'mse')

setDT(unimodal.sim)[, id := .GRP, by = .(seed)]

unimodal.sim$actual_start <- floor(unimodal.sim$peak_min)
unimodal.sim$actual_end <- ceiling(unimodal.sim$peak_max)
unimodal.sim <- unimodal.sim[!is.na(unimodal.sim$start), ]

jaccard.args <- unimodal.sim[, c('start', 'end', 'actual_start', 'actual_end')]

# Compute jaccard on the segment
unimodal.sim$overlap_size <- do.call(
  mapply,
  c(
    FUN = overlap_size,
    setNames(jaccard.args, c('a1', 'a2', 'b1', 'b2'))
    )
  )

unimodal.sim$union_size <- do.call(
  mapply,
  c(
    FUN = union_size,
    setNames(jaccard.args, c('a1', 'a2', 'b1', 'b2'))
    )
  )
rm(jaccard.args)

# Compute F(b) - F(a) to get a better representation of segment accuracy
unimodal.sim$prob_segment <- unlist(lapply(seq_along(unimodal.sim$actual_dist), function(i) {
  actual_dist <- unimodal.sim$actual_dist[[i]]
  param <- unimodal.sim$param[i]
  arg.params <- if (actual_dist == 'norm') {
    list(
      mean = 0,
      sd = param
    )
  } else if (actual_dist == 'gamma') {
    list(
      shape = param,
      rate = 1
    )
  } else if (actual_dist == 'unif') {
    list(
      min = unimodal.sim$actual_start[[i]],
      max = unimodal.sim$actual_end[[i]]
    )
  }

  segment_prob(
    distribution = actual_dist,
    params = arg.params,
    a = unimodal.sim$start[[i]],
    b = unimodal.sim$end[[i]]
    )
  }))

# unimodal.sim$seg_length == unimodal.sim$actu
unimodal.sim$actual_length <- unimodal.sim$actual_end - unimodal.sim$actual_start

unimodal.sim$jaccard <- unimodal.sim$overlap_size / unimodal.sim$union_size

unimodal.sim$correct_dist <- unimodal.sim$actual_dist == unimodal.sim$dist

# quantile_p <- seq(0, 1, by = 0.1) # c(0, 1/4, 2/4, 3/4, 1)
quantile_p <- c(0, 1/4, 2/4, 3/4, 1)
decile_p <- seq(0, 1, by = 0.1)

# Note: This are improperly called deciles even though they are quantiles
unimodal.sim$noise_decile <- cut(unimodal.sim$noise, scale.param(quantile_p, param.range = unimodal.params$noise))
unimodal.sim$eps_decile <- cut(unimodal.sim$eps, scale.param(quantile_p, param.range = unimodal.params$eps))
unimodal.sim$N_decile <- cut(unimodal.sim$N, round(scale.param(quantile_p, param.range = unimodal.params$N)))

# Only keep the segment with largest Jaccard
best.segment <- unimodal.sim[
  unimodal.sim[metric %in% c('mle', 'ks'), .I[which.max(jaccard)], by=.(id, metric)]$V1
  ]

unimodal.sim <- unimodal.sim[unimodal.sim[, .I[which.max(jaccard)], by=.(id, metric)]$V1]

unimodal.sim$jaccard_decile <- cut(unimodal.sim$jaccard, quantile(unimodal.sim$jaccard, p = decile_p))
unimodal.sim$actual_length_decile <- cut(unimodal.sim$actual_length, quantile(unimodal.sim$actual_length, p = decile_p))

unimodal.sim[
  ,
  .(seg_length = mean(actual_length),
    seg_length_median = median(actual_length),
    max_seg_length = max(actual_length),
    min_seg_length = min(seg_length),
    mean_param1 = mean(dist_param1, na.rm = TRUE),
    mean_param2 = mean(dist_param2, na.rm = TRUE),
    mean_true_param = mean(param, na.rm = TRUE)
    )
  , by = .(actual_dist)
  ]

### Create heatmap of the correlations
## TODO: Put this in a function
jaccard.cor.data <- as.data.frame(unimodal.sim[
  ,
  .(
    N_cor = cor(N, jaccard, method = 'spearman'),
    eps_cor = cor(eps, jaccard, method = 'spearman'),
    noise_cor = cor(noise, jaccard, method = 'spearman'),
    param_cor = cor(param, jaccard, method = 'spearman')
    )
  , by = .(remove_low_entropy, max_uniform, actual_dist)
  ])

cor.cols <- colnames(jaccard.cor.data)[grepl('_cor', colnames(jaccard.cor.data))]

jaccard.cor.order <- diana(jaccard.cor.data[, cor.cols])$order
jaccard.cor.order.t <- diana(t(jaccard.cor.data[, cor.cols]))$order

cor.cov.heatmap <- sim.plot.heatmap.cov(
    jaccard.cor.data[jaccard.cor.order, c(
      'remove_low_entropy', 'max_uniform', 'actual_dist'
      )]
    );

cor.data <- jaccard.cor.data[jaccard.cor.order, cor.cols[jaccard.cor.order.t]]

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
  )[jaccard.cor.order.t]

cov.hm <- create.heatmap(
  x = cor.data,
  same.as.matrix = TRUE,
  colourkey.cex = 2,
  clustering.method = 'none',
  row.pos = row.col.text$row,
  col.pos = row.col.text$col,
  cell.text = cor.data.round[abs(cor.data) > 0.1],
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
  main = 'Spearman correlation with segment Jaccard',
  main.cex = 2,
  layout.width = 2,
  layout.height = 1,
  resolution = 400,
  legend = list(
      right = list(
        fun = common.sim.legend(
          include.legends = c('params', 'distributions'),
          params.to.include = c('max_uniform', 'remove_low_entropy'),
          cont.params.to.include = NULL
          )
        )
      ),
  xlab.label = 'Spearman correlation',
  filename = print(
    file.path(
      plots.folder,
      generate.filename(paste0('HZSimulation', sim.version), 'cor-jaccard-segment', 'png')
      )
    )
  )

# cor_results_jaccard <- lapply(unimodal.sim[, c('N', 'noise', 'eps')], cor, method = 'spearman', y = unimodal.sim$jaccard)


## Evaulation of segments
## Jaccard overlap plots
sim.plot.segment.eval(
  best.segment,
  cluster = FALSE,
  print.colour.key = FALSE,
  xlab.label = 'Median Jaccard',
  resolution = 200,
  filename = print(
    file.path(
      plots.folder,
      generate.filename(paste0('HZSimulation', sim.version), 'median-jaccard', 'png')
      )
    )
  )

sim.plot.segment.eval(
  best.segment,
  cluster = TRUE,
  resolution = 200,
  xlab.label = 'Median Jaccard',
  filename = print(
    file.path(
      plots.folder,
      generate.filename(paste0('HZSimulation', sim.version), 'median-jaccard-clustered', 'png')
      )
    )
  )

sim.plot.segment.eval(
  best.segment,
  cluster = FALSE,
  target = 'prob_segment',
#   xlab.label = 'Median probability of segment',
  print.colour.key = TRUE,
  resolution = 200,
  filename = print(
    file.path(
      plots.folder,
      generate.filename(paste0('HZSimulation', sim.version), 'median-CDF-prob', 'png')
      )
    )
  )

sim.plot.segment.eval(
  best.segment,
  cluster = TRUE,
  target = 'prob_segment',
  xlab.label = 'Median probability of segment',
  print.colour.key = TRUE,
  resolution = 200,
  filename = print(
    file.path(
      plots.folder,
      generate.filename(paste0('HZSimulation', sim.version), 'median-CDF-prob-clustered', 'png')
      )
    )
  )

# Create the probability of segment histogram plot
prob.segment.hists.rle <- lapply(c('gamma', 'norm', 'unif'), function(d) {
  dist_name <- switch (d,
    'unif' = 'Uniform',
    'gamma' = 'Gamma',
    'norm' = 'Normal'
    )
  lapply(c(FALSE, TRUE), function(remove_low_entropy) {
    breaks <- seq(0, 1, by = 0.05)
    yat <- seq(0, 100, by = 20)
    create.histogram(
      best.segment$prob_segment[best.segment$actual_dist == d & best.segment$remove_low_entropy == remove_low_entropy],
      breaks = breaks,
      xlimits = c(0, 1),
      ylimits = c(0, 110),
      yat = yat,
      col = adjustcolor(distribution_colours[[d]], alpha.f = 0.4),
      ylab.label = if (!remove_low_entropy) paste0(dist_name,'\nPercent') else '\n',
      yaxis.lab = if (!remove_low_entropy) TRUE else rep(' ', length(yat)),
      xlab.label = if (d == 'unif') as.character(remove_low_entropy) else '',
      yaxis.tck = if (!remove_low_entropy) c(1, 0) else 0,
      xaxis.tck = if (d == 'unif') c(1, 0) else 0,
      xaxis.lab = if (d == 'unif')  TRUE else rep('', length(breaks)),
    )
  })
})

create.multipanelplot(
  unlist(prob.segment.hists.rle, recursive = FALSE),
  layout.height = length(prob.segment.hists.rle),
  layout.width = 2,
  main = 'Probability of segment',
  xlab.label = 'Remove low entropy',
  xlab.cex = 2.5,
  x.spacing = -1.5,
  main.cex = 2.5,
  width = 12,
  height = 12,
  resolution = 300,
  filename = print(
    file.path(
      plots.folder,
      generate.filename(paste0('HZSimulation', sim.version), 'prob-segment-histogram', 'png')
      )
    )
  )

## Distribution accuracy

acc <- 'dist'

categorical.vars <- c(
      'actual_dist', 'max_uniform', 'remove_low_entropy'
    )
decile.vars <- c('eps_decile', 'N_decile', 'noise_decile', 'jaccard_decile')

decile.uni.plots <- lapply(decile.vars, function(v) {
  main <- sub('_decile', '', v)
  res <- sim.plot.quantile.accuracy(
    unimodal.sim,
    resolution = 200,
    cluster = FALSE,
    group_vars = c(categorical.vars, v),
    acc = acc,
    legend = NULL,
    main = main,
    print.colour.key = FALSE,
    xlab.label = ''
    )
    class(res) <- c('frame', 'gTree', 'grob', 'gDesc')
    res
  })

decile.uni.plots.clustered <- lapply(decile.vars, function(v) {
  main <- sub('_decile', '', v)
  res <- sim.plot.quantile.accuracy(
    unimodal.sim,
    resolution = 200,
    cluster = TRUE,
    group_vars = c(categorical.vars, v),
    acc = acc,
    legend = NULL,
    main = main,
    print.colour.key = FALSE,
    xlab.label = ''
    )
    class(res) <- c('frame', 'gTree', 'grob', 'gDesc')
    res
  })

decile.uni.plots.clustered <- lapply(decile.vars, function(v) {
  main <- sub('_decile', '', v)
  res <- sim.plot.quantile.accuracy(
    unimodal.sim,
    resolution = 200,
    cluster = TRUE,
    group_vars = c(categorical.vars, v),
    acc = acc,
    legend = NULL,
    main = main,
    print.colour.key = FALSE,
    xlab.label = ''
    )
    class(res) <- c('frame', 'gTree', 'grob', 'gDesc')
    res
  })

# This is bug in BPG
# Need to make a fake device for colourkey to write to
png(filename = tempfile(fileext = '.png'))
colourkey <- create.colourkey(
    x = seq(0, 1, length.out = 100),
    at = seq(0, 1, length.out = 20),
    colour.scheme = c('white', 'red'),
    colourkey.labels.cex = 2.5,
    placement = viewport(
        just='center',
        x = 0.5,
        y = 1,
        width = 0.5,
        height = 0.1
        )
    );
dev.off()

create.multipanelplot(
  decile.uni.plots,
  width = 20,
  height = 20,
  layout.height = 2,
  layout.width = 2,
  resolution = 200,
  main = 'No Clustering',
  xlab.label = 'Distribution Accuracy',
  legend = list(
    right = list(
      fun = common.sim.legend(
        include.legends = c('params', 'distributions', 'quantiles'),
        params.to.include = c('max_uniform', 'remove_low_entropy'),
        cont.params.to.include = decile.vars
        )
      ),
    bottom = list(
      fun = colourkey
    )
  ),
  filename = print(
    file.path(
      plots.folder,
      generate.filename(paste0('HZSimulation', sim.version), 'dist-acc-unimodal-mpp', 'png')
      )
    )
  )

create.multipanelplot(
  decile.uni.plots.clustered,
  width = 20,
  height = 20,
  layout.height = 2,
  layout.width = 2,
  resolution = 200,
  main = 'DIANA Clustered',
  xlab.label = 'Distribution Accuracy',
  legend = list(
    right = list(
      fun = common.sim.legend(
        include.legends = c('params', 'distributions', 'quantiles'),
        params.to.include = c('max_uniform', 'remove_low_entropy'),
        cont.params.to.include = decile.vars
        )
      ),
    bottom = list(
      fun = colourkey
    )
  ),
  filename = print(
    file.path(
      plots.folder,
      generate.filename(paste0('HZSimulation', sim.version), 'dist-acc-unimodal-mpp-cluster', 'png')
      )
    )
  )

# Across all variables
sim.plot.quantile.accuracy(
    unimodal.sim,
    resolution = 800,
    cluster = TRUE,
    acc = 'dist',
    group_vars = c(
      'actual_dist', 'max_uniform', 'remove_low_entropy',
      'N_decile', 'noise_decile'
      ),
    main = 'DIANA Clustered Distribution Accuracy',
    print.colour.key = TRUE,
    xlab.label = 'Distribution Accuracy',
    filename = print(
    file.path(
      plots.folder,
      generate.filename(paste0('HZSimulation', sim.version), 'dist-acc-unimodal-mpp-cluster-all', 'png')
      )
    )
  )

sim.plot.quantile.accuracy(
    unimodal.sim,
    resolution = 400,
    cluster = FALSE,
    acc = 'dist',
    main = 'No clustering',
    group_vars = c(
      'actual_dist', 'max_uniform', 'remove_low_entropy',
      'N_decile', 'noise_decile'
      ),
    print.colour.key = TRUE,
    xlab.label = 'Distribution Accuracy',
    filename = print(
    file.path(
      plots.folder,
      generate.filename(paste0('HZSimulation', sim.version), 'dist-acc-unimodal-mpp-all', 'png')
      )
    )
  )

### Prints the alternative metrics on the Noise and N quantiles
for (eval_metric in c('F1', 'Precision', 'Recall', 'Sensitivity', 'Specificity', 'Balanced Accuracy')) {
  no.clust.cm <- sim.plot.quantile.accuracy(
    unimodal.sim,
    cluster = FALSE,
    acc = eval_metric,
    legend = NULL,
    main = paste0(eval_metric, ': No clustering'),
    group_vars = c(
        'max_uniform', 'remove_low_entropy',
        'N_decile', 'noise_decile'
        ),
      print.colour.key = FALSE,
      xlab.label = ''
  )
  class(no.clust.cm) <- c('multipanel', 'frame', 'gTree', 'grob', 'gDesc')

  clust.cm <- sim.plot.quantile.accuracy(
    unimodal.sim,
    cluster = TRUE,
    acc = eval_metric,
    legend = NULL,
    main = paste0(eval_metric, ': DIANA Clustered'),
    group_vars = c(
        'max_uniform', 'remove_low_entropy',
        'N_decile', 'noise_decile'
        ),
      print.colour.key = FALSE,
      xlab.label = ''
    )
  class(clust.cm) <- c('multipanel', 'frame', 'gTree', 'grob', 'gDesc')

  create.multipanelplot(
  list(no.clust.cm, clust.cm),
    width = 20,
    height = 10,
    layout.height = 1,
    layout.width = 2,
    resolution = 500,
    main = '',
    xlab.label = eval_metric,
    legend = list(
      right = list(
        fun = common.sim.legend(
          include.legends = c('params', 'distributions', 'quantiles'),
          params.to.include = c('max_uniform', 'remove_low_entropy'),
          cont.params.to.include = c('N_decile', 'noise_decile')
          )
        ),
      bottom = list(
        fun = colourkey
      )
    ),
    filename = print(
      file.path(
        plots.folder,
        generate.filename(paste0('HZSimulation', sim.version), paste0(sub(' ', '-', eval_metric), '-unimodal-mpp-all'), 'png')
        )
      )
  )
}

