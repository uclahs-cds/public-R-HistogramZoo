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
  pattern = '2023-07-18_.*v5.*noise.tsv',
  full.names = TRUE
  )

unimodal.sim.metrics <- rbindlist(
  lapply(metric.files, data.table::fread, sep = '\t'),
  fill = TRUE
  )

mle.files <- list.files(
  path = merged.folder,
  pattern = '2023-07-18_.*v5.*mle',
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
# unimodal.sim$seg_length == unimodal.sim$actu
unimodal.sim$actual_length <- unimodal.sim$actual_end - unimodal.sim$actual_start

unimodal.sim$jaccard <- unimodal.sim$overlap_size / unimodal.sim$union_size

# hist(unimodal.sim$jaccard)
# Remove really high noise
# unimodal.sim <- unimodal.sim[noise <= 0.5, ]

# Compute Jaccard

# (unimodal.sim[abs(unimodal.sim$noise_min - unimodal.sim$noise_min.validation) > 1, ])
#
# mean(abs(unimodal.sim$noise_min - unimodal.sim$noise_min.validation) > 1)
#
# View(unimodal.sim[abs(unimodal.sim$noise_min - unimodal.sim$noise_min.validation) > 1, ])

## Segment level analysis:
# Jaccard (compare with actual segment)

# quantile_p <- seq(0, 1, by = 0.1) # c(0, 1/4, 2/4, 3/4, 1)
quantile_p <- c(0, 1/4, 2/4, 3/4, 1)

unimodal.sim$jaccard_decile <- cut(unimodal.sim$jaccard, breaks = seq(0, 1, by = 0.1))
unimodal.sim$noise_decile <- cut(unimodal.sim$noise, scale.param(quantile_p, param.range = unimodal.params$noise))
unimodal.sim$eps_decile <- cut(unimodal.sim$eps, scale.param(quantile_p, param.range = unimodal.params$eps))
unimodal.sim$N_decile <- cut(unimodal.sim$N, round(scale.param(quantile_p, param.range = unimodal.params$N)))

# Only keep the segment with largest Jaccard
best.segment <- unimodal.sim[
  unimodal.sim[, .I[which.max(jaccard)], by=.(id, metric)]$V1
  ]

# Which metric we select is arbitrary for union, ks, jaccard, etc since they have the same segmentation
best.segment <- best.segment[best.segment$metric %in% c('mle', 'ks'), ];

best.segment[
  ,
  .(jaccard = mean(jaccard),
    N = .N)
  ,by = .(dist, actual_dist),]


unimodal.sim <- unimodal.sim[unimodal.sim[, .I[which.max(jaccard)], by=.(id, metric)]$V1]

# View(unimodal.sim[, c('jaccard', 'start', 'end', 'actual_start', 'actual_end', 'union_size')])

sim.plot.segment.jaccard(
  best.segment,
  cluster = FALSE,
  resolution = 200,
  filename = print(
    file.path(
      plots.folder,
      generate.filename('HZSimulation', 'median-jaccard', 'png')
      )
    )
  )

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

for (acc in c('dist', 'peaks')) {
# for (acc in c('count')) {
  sim.plot.quantile.accuracy(
    unimodal.sim,
    resolution = 200,
    sort.cols = 'jaccard_decile',
    acc = acc,
    filename =  print(
      file.path(
        plots.folder,
        generate.filename('HZSimulation', paste0('quartile-classification-heatmap-', acc), 'png')
        )
      )
    )
}

# Sort

# Decile partitioned by parameters
for (m in metrics) {
  sim.plot.quantile.accuracy(
  unimodal.sim,
  acc = 'both',
  sort.cols = m,
  resolution = 200,
  filename =  print(
    file.path(
      plots.folder,
      generate.filename(
        'HZSimulation',
        paste0('decile-sort-classification-heatmap-', m),
        'png'
        )
      )
    )
  )
}

