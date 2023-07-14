library(HistogramZoo);
library(data.table);
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

base.path <- HistogramZoo:::load.config()$root.path;
results.folder <- file.path(base.path, 'results');
merged.folder <- file.path(results.folder, 'merged_sims');
plots.folder <- file.path(base.path, 'plots');

#sim.tsv.mle.paths <- list.files(results.folder, pattern = 'Unimodal_Sim_MLE.tsv$', full.names = TRUE);
#sim.tsv.paths <- list.files(results.folder, pattern = 'Unimodal_Sim.tsv$', full.names = TRUE);

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
# unimodal.sim <- unimodal.sim.metrics;
unimodal.sim <- rbindlist(
  list(
    unimodal.sim.metrics,
    unimodal.sim.metrics.mle
    ),
  fill = TRUE
  )

metrics <- c("mle", "chisq", "intersection", "jaccard", "ks", "mse")

# setDT(unimodal.sim)[, id := .GRP, by = .(N, param, noise, actual_dist, eps, metric)]
setDT(unimodal.sim)[, id := .GRP, by = .(seed)]

# Remove really high noise
unimodal.sim <- unimodal.sim[noise <= 0.5, ]

unimodal.sim$seg_length <-  unimodal.sim$end - unimodal.sim$start;

# Only keep the largest segment
unimodal.sim <- unimodal.sim[unimodal.sim[, .I[which.max(seg_length)], by=.(id, metric)]$V1]

# sim.data.largest$num_segments

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
unimodal.sim$noise_decile <- cut(unimodal.sim$noise, quantile(unimodal.sim$noise, probs = quantile_p))
unimodal.sim$eps_decile <- cut(unimodal.sim$eps, quantile(unimodal.sim$eps, probs = quantile_p))
unimodal.sim$N_decile <- cut(unimodal.sim$N, quantile(unimodal.sim$N, probs = quantile_p))

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
