library(HistogramZoo);
library(data.table);
library(BoutrosLab.utilities);

base.path <- HistogramZoo:::load.config()$root.path;
# sim.file <- file.path(base.path, 'results', 'unimodal_sim_noise', 'processed_results.tsv');
results.folder <- file.path(base.path, 'results');
merged.folder <- file.path(results.folder, 'merged_sims');

# v1-v3 need to have peak_min, peak_max computed
# sim.folders <- c(
#     'unimodal_sim_noise_v2', 'unimodal_sim_noise_v3',
#     'unimodal_sim_noise_v4', 'unimodal_sim_noise_v5');

sim.folders <- 'unimodal_sim_noise_v5';

for (sim in sim.folders) {
  for (mle in c(FALSE, TRUE)) {
    cat('Processing: ', sim, ' files... MLE = ', mle, '\n');
    sim.folder <- file.path(results.folder, sim);
    mle_pattern <- if (mle) 'Unimodal_Sim_MLE.tsv$' else 'Unimodal_Sim.tsv$'

    sim.tsv.paths <- list.files(sim.folder, pattern = mle_pattern, full.names = TRUE, include.dirs = FALSE);

    sim.data <- rbindlist(
        lapply(sim.tsv.paths, data.table::fread, sep = '\t'),
        fill = TRUE
        )

    sim.data$seg_length <-  sim.data$end - sim.data$start;

    print(head(sim.data))

    # Errors have newlines that can complicate things
    sim.data$error <- gsub('\n', '', sim.data$error);
    mle.suffix <- if (mle) '-mle' else '';

    write.table(
        sim.data,
        file = print(file.path(
            merged.folder,
            generate.filename(
                'HZSimulation',
                paste0(sim, '-merged-unimodal-sim-noise', mle.suffix),
                'tsv'
                )
            )),
        sep = '\t',
        row.names = FALSE
        )
    }
}


