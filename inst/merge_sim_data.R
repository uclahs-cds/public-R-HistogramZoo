library(HistogramZoo);
library(data.table);
library(BoutrosLab.utilities);

base.path <- HistogramZoo:::load.config()$root.path;
# sim.file <- file.path(base.path, 'results', 'unimodal_sim_noise', 'processed_results.tsv');

sim.folders <- c('unimodal_sim_noise_v2');

for (sim in sim.folders) {
    cat('Processing: ', sim, ' files...\n');

    results.folder <- file.path(base.path, 'results');
    sim.folder <- file.path(results.folder, sim);
    merged.folder <- file.path(results.folder, 'merged_sims');
    plots.folder <- file.path(base.path, 'plots');

    sim.tsv.mle.paths <- list.files(sim.folder, pattern = 'Unimodal_Sim_MLE.tsv$', full.names = TRUE);
    # sim.tsv.mle.paths <- sim.tsv.mle.paths[c(1, length(sim.tsv.mle.paths))];
    sim.tsv.paths <- list.files(sim.folder, pattern = 'Unimodal_Sim.tsv$', full.names = TRUE);

    sim.mle.data <- rbindlist(
        lapply(sim.tsv.mle.paths, data.table::fread, sep = '\t'),
        fill = TRUE
        )

    sim.data <- rbindlist(
        lapply(sim.tsv.paths, data.table::fread, sep = '\t'),
        fill = TRUE
        )

    sim.data$seg_length <-  sim.data$end - sim.data$start;

    write.table(
        sim.data,
        file = print(file.path(
            merged.folder,
            generate.filename(
                'HZSimulation',
                'merged-unimodal-sim-noise',
                'tsv'
                )
            )),
        sep = '\t',
        row.names = FALSE
        )

    write.table(
        sim.mle.data,
        file = print(file.path(
            merged.folder,
            generate.filename(
                'HZSimulation',
                'merged-unimodal-sim-noise-MLE',
                'tsv'
                )
            )),
        sep = '\t',
        row.names = FALSE
        )
}


