library(HistogramZoo);
library(data.table);

base.path <- HistogramZoo:::load.config()$root.path;
sim.file <- file.path(base.path, 'results', 'unimodal_sim_noise', 'processed_results.tsv');
# tsv.paths <- list.files(results.folder, pattern = '_Unimodal_Sim[.]rds_processed_results[.]tsv', full.names = TRUE);
sim.data <- fread(sim.file, sep = '\t')

# Determine accuracy
# print(read.table(tsv.paths[[1]], sep = '\t'))
dim(sim.data)


sim.data$seg_length <-  sim.data$end - sim.data$start
# Only keep the largest segment
sim.data <- sim.data[sim.data[, .I[which.max(seg_length)], by=.(id, metric)]$V1]

sim.data[, .(correct = actual_dist == dist), by = .(id, metric)
    ][
    , .(accuracy = mean(correct, na.rm = TRUE)), by = .(metric)
    ]

