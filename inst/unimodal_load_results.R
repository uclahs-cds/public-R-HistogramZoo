library(HistogramZoo);
library(data.table);

base.path <- HistogramZoo:::load.config()$root.path;
results.folder <- file.path(base.path, 'results', 'unimodal_sim_noise');
sim.tsv.paths <- list.files(results.folder, pattern = 'Unimodal_Sim.tsv$', full.names = TRUE);

all_results <- rbindlist(
    lapply(sim.tsv.paths, data.table::fread, sep = '\t'),
    fill = TRUE
    )

filename <- file.path(results.folder, 'processed_results.tsv');
cat('Writing to: ', filename, '\n');

# TODO: Split this into batches and combine?
write.table(
    x = all_results,
    file = file.path(results.folder, 'processed_results.tsv'),
    sep = '\t',
    row.names = FALSE
    )
