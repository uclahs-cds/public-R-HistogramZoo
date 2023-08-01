library(HistogramZoo);
library(data.table);
library(BoutrosLab.utilities);

base.path <- HistogramZoo:::load.config()$root.path;
results.folder <- file.path(base.path, 'results');
merged.folder <- file.path(results.folder, 'merged_sims');

sim.folders <- paste0('multi_modal_sim_v', 3);

sim.file.type <- c('overlap', 'seg_results', 'actual_peaks');

# It might be better to do individual processing (merging) before combining
# For now it just merges all sim runs from each file type into combined file
for (sim in sim.folders) {
  for (file.type in sim.file.type) {
    for (mle in c(FALSE, TRUE)) {
      cat('Processing: ', sim, ' ', file.type, ' files... MLE = ', mle, '\n');
      sim.folder <- file.path(results.folder, sim);

      file.pattern <- '[0-9]+_multi-peak-sim_'
      if (mle)  file.pattern <- paste0(file.pattern, '_?MLE_')
      file.pattern <- paste0(file.pattern, file.type, '[.]tsv')

      sim.tsv.paths <- list.files(sim.folder, pattern = file.pattern, full.names = TRUE, include.dirs = FALSE);

      sim.data <- rbindlist(
          lapply(sim.tsv.paths, data.table::fread, sep = '\t'),
          fill = TRUE
          )

      if (nrow(sim.data) == 0) next;

      suffix <- if (mle) '-mle' else '';

      write.table(
          sim.data,
          file = print(file.path(
              merged.folder,
              generate.filename(
                  'HZSimulation',
                  paste0(sim, '-merged-multi-modal-', file.type, suffix),
                  'tsv'
                  )
              )),
          sep = '\t',
          row.names = FALSE
          )
      }
    }
  }


