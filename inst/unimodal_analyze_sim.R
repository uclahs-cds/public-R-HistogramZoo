library(HistogramZoo);
library(data.table);
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

base.path <- HistogramZoo:::load.config()$root.path;
results.folder <- file.path(base.path, 'results');
merged.folder <- file.path(results.folder, 'merged_sims');
plots.folder <- file.path(base.path, 'plots');


unimodal.sim.metrics <- data.table::fread(
  file = file.path(
    merged.folder,
    '2023-07-13_HZSimulation_merged-unimodal-sim-noise.tsv'
    ),
  sep = '\t'
  )

unimodal.sim.metrics.mle <- data.table::fread(
  file = file.path(
    merged.folder,
    '2023-07-13_HZSimulation_merged-unimodal-sim-noise-MLE.tsv'
    ),
  sep = '\t'
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

# Remove really high noise
unimodal.sim <- unimodal.sim[noise <= 0.5, ]

unimodal.sim$seg_length <-  unimodal.sim$end - unimodal.sim$start;

# Only keep the largest segment
unimodal.sim <- unimodal.sim[unimodal.sim[, .I[which.max(seg_length)], by=.(id, metric)]$V1]

# TODO: Need to generate the data with the given seed, then can get the peak_min, peak_max
# Note: This is temporary until the new simulations

for (acc in c('dist', 'peaks', 'both')) {
  sim.plot.overall.accuracy(
    unimodal.sim,
    resolution = 200,
    acc = acc,
    filename = print(
      file.path(
        plots.folder,
        generate.filename(
          'HZSimulation',
          paste0('overall-classification-heatmap-', acc),
          'png'
          )
        )
      )
    )
}

quantile_p <- seq(0, 1, by = 0.1) # c(0, 1/4, 2/4, 3/4, 1)

scale.param <- function(x, param.range) {
  (param.range[2] - param.range[1]) * x + param.range[1]
}

unimodal.sim$noise_decile <- cut(unimodal.sim$noise, scale.param(quantile_p, param.range = unimodal.params$noise))
unimodal.sim$eps_decile <- cut(unimodal.sim$eps, scale.param(quantile_p, param.range = unimodal.params$eps))
unimodal.sim$N_decile <- cut(unimodal.sim$N, round(scale.param(quantile_p, param.range = unimodal.params$N)))

for (acc in c('dist', 'peaks')) {
  sim.plot.quantile.accuracy(
    unimodal.sim,
    resolution = 200,
    acc = acc,
    filename =  print(
      file.path(
        plots.folder,
        generate.filename('HZSimulation', paste0('decile-classification-heatmap-', acc), 'png')
        )
      )
    )
  }

# Decile partitioned by parameters
for (m in metrics) {
  sim.plot.quantile.accuracy(
  unimodal.sim,
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

