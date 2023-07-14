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

setDT(unimodal.sim)[, id := .GRP, by = .(N, param, noise, actual_dist, eps, metric)]

# Remove really high noise
unimodal.sim <- unimodal.sim[noise <= 0.5, ]

unimodal.sim$seg_length <-  unimodal.sim$end - unimodal.sim$start;

# Only keep the largest segment
unimodal.sim <- unimodal.sim[unimodal.sim[, .I[which.max(seg_length)], by=.(id, metric)]$V1]

sim.plot.overall.accuracy(
  unimodal.sim,
  resolution = 200,
  filename = print(
    file.path(
      plots.folder,
      generate.filename('HZSimulation', 'overall-classification-heatmap', 'png')
      )
    )
  );
