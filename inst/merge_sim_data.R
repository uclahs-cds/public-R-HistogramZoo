library(HistogramZoo);
library(data.table);
library(BoutrosLab.utilities);

base.path <- HistogramZoo:::load.config()$root.path;
# sim.file <- file.path(base.path, 'results', 'unimodal_sim_noise', 'processed_results.tsv');
results.folder <- file.path(base.path, 'results');
merged.folder <- file.path(results.folder, 'merged_sims');

# v1-v3 need to have peak_min, peak_max computed
sim.folders <- c('unimodal_sim_noise', 'unimodal_sim_noise_v2', 'unimodal_sim_noise_v3');

for (sim in sim.folders) {
  for (mle in c(FALSE, TRUE)) {
    cat('Processing: ', sim, ' files... MLE = ', mle, '\n');
    sim.folder <- file.path(results.folder, sim);
    mle_pattern <- if (mle) 'Unimodal_Sim_MLE.tsv$' else 'Unimodal_Sim.tsv$'

    sim.tsv.paths <- list.files(sim.folder, pattern = , full.names = TRUE);

    sim.data <- rbindlist(
        lapply(sim.tsv.paths, data.table::fread, sep = '\t'),
        fill = TRUE
        )

    sim.data$seg_length <-  sim.data$end - sim.data$start;

    if (! c('peak_min') %in% colnames(sim.data)) {
      unique.recompute <- unique(
      sim.data[, c('seed', 'actual_dist', 'max_uniform', 'remove_low_entropy', 'truncated_models')]
      )

    recompute.peak.rage <- do.call(
      'rbind.data.frame',
      mapply(
        HistogramZoo:::peak_min_recompute,
        seed = unique.recompute$seed,
        actual_dist = unique.recompute$actual_dist,
        max_uniform = unique.recompute$max_uniform,
        remove_low_entropy = unique.recompute$remove_low_entropy,
        SIMPLIFY = FALSE
        )
      )

    sim.data <- merge(
      x = sim.data,
      y = recompute.peak.rage,
      by = 'seed',
      suffixes = c('', '.validation')
      )
    }

    mle.suffix <- if (mle) '-mle' else '';

    write.table(
        sim.data,
        file = print(file.path(
            merged.folder,
            generate.filename(
                'HZSimulation',
                paste0('merged-unimodal-sim-noise', mle.suffix),
                'tsv'
                )
            )),
        sep = '\t',
        row.names = FALSE
        )
    }
}


