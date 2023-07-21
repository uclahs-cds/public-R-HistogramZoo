library(HistogramZoo);
library(data.table);
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

base.path <- HistogramZoo:::load.config()$root.path;
results.folder <- file.path(base.path, 'results');
merged.folder <- file.path(results.folder, 'merged_sims');
plots.folder <- file.path(base.path, 'plots');

metric.files <- list.files(
  path = merged.folder,
  pattern = '2023-07-21_.*v5.*noise.tsv',
  full.names = TRUE
  )

unimodal.sim.metrics <- rbindlist(
  lapply(metric.files, data.table::fread, sep = '\t'),
  fill = TRUE
  )

mle.files <- list.files(
  path = merged.folder,
  pattern = '2023-07-21_.*v5.*mle',
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
# unimodal.sim$seg_length == unimodal.sim$actu
unimodal.sim$actual_length <- unimodal.sim$actual_end - unimodal.sim$actual_start

unimodal.sim$jaccard <- unimodal.sim$overlap_size / unimodal.sim$union_size

unimodal.sim$correct_dist <- unimodal.sim$actual_dist == unimodal.sim$dist

# quantile_p <- seq(0, 1, by = 0.1) # c(0, 1/4, 2/4, 3/4, 1)
quantile_p <- c(0, 1/4, 2/4, 3/4, 1)
decile_p <- seq(0, 1, by = 0.1)

unimodal.sim$jaccard_decile <- cut(unimodal.sim$jaccard, quantile(unimodal.sim$jaccard, p = decile_p))
unimodal.sim$noise_decile <- cut(unimodal.sim$noise, scale.param(quantile_p, param.range = unimodal.params$noise))
unimodal.sim$eps_decile <- cut(unimodal.sim$eps, scale.param(quantile_p, param.range = unimodal.params$eps))
unimodal.sim$N_decile <- cut(unimodal.sim$N, round(scale.param(quantile_p, param.range = unimodal.params$N)))

# Only keep the segment with largest Jaccard
best.segment <- unimodal.sim[
  unimodal.sim[metric %in% c('mle', 'ks'), .I[which.max(jaccard)], by=.(id, metric)]$V1
  ]

unimodal.sim <- unimodal.sim[unimodal.sim[, .I[which.max(jaccard)], by=.(id, metric)]$V1]

unimodal.sim[
  ,
  .(seg_length = mean(actual_length),
    seg_length_median = median(actual_length),
    max_seg_length = max(actual_length),
    min_seg_length = min(seg_length)),
  , by = .(actual_dist)
  ]

# Which metric we select is arbitrary for union, ks, jaccard, etc since they have the same segmentation
# best.segment <- best.segment[best.segment$metric %in% c('mle', 'ks'), ];

# class.acc <- unimodal.sim[
#   dist == actual_dist,
#   .(jaccard = mean(jaccard),
#     N = .N)
#   , by = .(dist, actual_dist, metric),
#   ]
# class.acc$perc <- class.acc$N / sum(class.acc$N)


# dist.log.mod <- glm(
#   formula = correct_dist ~ actual_dist + N + actual_dist:param + eps +
#     noise + max_uniform + remove_low_entropy + metric + actual_dist:jaccard,
#   data = unimodal.sim,
#   family = binomial()
#   )
# dist.log.mod.coefs <- as.data.frame(summary(dist.log.mod)$coefficients)
# dist.log.mod.coefs$est.exp <- exp(dist.log.mod.coefs$Estimate)

# kableExtra::kable_styling(kableExtra::kable(dist.log.mod.coefs))

# View(unimodal.sim[, c('jaccard', 'start', 'end', 'actual_start', 'actual_end', 'union_size')])

# create.scatterplot(
#     formula = jaccard ~ eps | actual_dist,
#     data = best.segment,
#     resolution = 200,
#     alpha = 0.2,
#     filename = 'test_jaccard_eps_scatter.png'
#   )
#
# cor(best.segment$eps, best.segment$jaccard, method = 'spearman')
# cor(best.segment$eps, best.segment$num_segments, method = 'spearman')
#
# table(best.segment$num_segments)
# cor(jaccard ~ eps, best.segment)

sim.plot.segment.jaccard(
  best.segment,
  cluster = FALSE,
  print.colour.key = FALSE,
  resolution = 200,
  filename = print(
    file.path(
      plots.folder,
      generate.filename('HZSimulation', 'median-jaccard', 'png')
      )
    )
  )

# Shouldn't really have major differences since we are sampling uniformly
# sim.plot.segment.jaccard(
#   best.segment,
#   cluster = FALSE,
#   print.colour.key = FALSE,
#   resolution = 200,
#   target = 'count',
#   filename = print(
#     file.path(
#       plots.folder,
#       generate.filename('HZSimulation', 'count-jaccard', 'png')
#       )
#     )
#   )

sim.plot.segment.jaccard(
  best.segment,
  cluster = TRUE,
  resolution = 200,
  filename = print(
    file.path(
      plots.folder,
      generate.filename('HZSimulation', 'median-jaccard-clustered', 'png')
      )
    )
  )

#for (cluster in c(FALSE, TRUE)){
# cluster <- FALSE
#for (acc in c('count', 'dist', 'peaks')) {
acc <- 'dist'

categorical_vars <- c(
      'actual_dist', 'max_uniform', 'remove_low_entropy'
    )
decile_vars <- c('eps_decile', 'N_decile', 'noise_decile', 'jaccard_decile')

decile_uni_plots <- lapply(decile_vars, function(v) {
  main <- sub('_decile', '', v)
  res <- sim.plot.quantile.accuracy(
    unimodal.sim,
    resolution = 200,
    cluster = FALSE,
    group_vars = c(categorical_vars, v),
    acc = acc,
    legend = NULL,
    main = main,
    print.colour.key = FALSE,
    xlab.label = ''
    )
    class(res) <- c('frame', 'gTree', 'grob', 'gDesc')
    res
  })

decile_uni_plots_clustered <- lapply(decile_vars, function(v) {
  main <- sub('_decile', '', v)
  res <- sim.plot.quantile.accuracy(
    unimodal.sim,
    resolution = 200,
    cluster = TRUE,
    group_vars = c(categorical_vars, v),
    acc = acc,
    legend = NULL,
    main = main,
    print.colour.key = FALSE,
    xlab.label = ''
    )
    class(res) <- c('frame', 'gTree', 'grob', 'gDesc')
    res
  })

png(nullfile())
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

create.multipanelplot(
  decile_uni_plots,
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
        cont.params.to.include = decile_vars
        )
      ),
    bottom = list(
      fun = colourkey
    )
  ),
  filename = print(
    file.path(
      plots.folder,
      generate.filename('HZSimulation', 'dist-acc-unimodal-mpp', 'png')
      )
    )
  )

create.multipanelplot(
  decile_uni_plots_clustered,
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
        cont.params.to.include = decile_vars
        )
      ),
    bottom = list(
      fun = colourkey
    )
  ),
  filename = print(
    file.path(
      plots.folder,
      generate.filename('HZSimulation', 'dist-acc-unimodal-mpp-cluster', 'png')
      )
    )
  )

# Scatter plot of # of segments vs eps?
#

