library(HistogramZoo);
library(data.table);
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

base.path <- HistogramZoo:::load.config()$root.path;
# sim.file <- file.path(base.path, 'results', 'unimodal_sim_noise', 'processed_results.tsv');
results.folder <- file.path(base.path, 'results', 'unimodal_sim_noise_mle');

sim.tsv.paths <- list.files(results.folder, pattern = 'Unimodal_Sim_MLE.tsv$', full.names = TRUE);

sim.data <- rbindlist(
    lapply(sim.tsv.paths, data.table::fread, sep = '\t'),
    fill = TRUE
    )


sim.data$seg_length <-  sim.data$end - sim.data$start

# write.table(sim.data, file = file.path(results.folder, 'unimodal_sim_noise_mle.tsv'), sep = '\t', row.names = FALSE)

if (FALSE) sim.data <- read.table(
  file = file.path(results.folder, 'unimodal_sim_noise_mle.tsv'),
  sep = '\t',
  header = TRUE
  )

setDT(sim.data)[, id := .GRP, by = .(N, param, noise, actual_dist, eps)]

# Only keep the largest segment
sim.data.largest <- sim.data[sim.data[, .I[which.max(seg_length)], by=.(id, metric)]$V1]

# sim.data.largest$num_segments

# Accuracy by distribution
sim.data.largest[
    , .(
      accuracy_dist = mean(actual_dist == dist, na.rm = TRUE),
      accuracy_peaks = mean(num_segments == 1, na.rm = TRUE),
      N = mean(N),
      eps = mean(eps),
      noise = mean(noise)
      ), by = .(metric, actual_dist)
    ]

# TODO: Deciles
quantile_p <- c(0, 1/4, 2/4, 3/4, 1)
sim.data.largest$noise_quartile <- cut(sim.data.largest$noise, quantile(sim.data.largest$noise, probs = quantile_p))
sim.data.largest$eps_quartile <- cut(sim.data.largest$eps, quantile(sim.data.largest$eps))
sim.data.largest$N_quartile <- cut(sim.data.largest$N, quantile(sim.data.largest$N))
levels(sim.data.largest$noise_quartile) <- quantile_p[-1]
levels(sim.data.largest$eps_quartile) <- quantile_p[-1]
levels(sim.data.largest$N_quartile) <- quantile_p[-1]

sim.data.largest$correct <- as.factor(sim.data.largest$actual_dist == sim.data.largest$dist)

sim.data.largest[
    order(
      actual_dist,
      N_quartile,
      noise_quartile,
      eps_quartile
    ), .(
      accuracy_dist = mean(actual_dist == dist, na.rm = TRUE),
      accuracy_peaks = mean(num_segments == 1, na.rm = TRUE),
      N = mean(N),
      eps = mean(eps),
      noise = mean(noise)
      ), by = .(metric, actual_dist, noise_quartile, eps_quartile, N_quartile)
    ]

sim.data.largest[
    , .(
      accuracy_dist = mean(actual_dist == dist, na.rm = TRUE),
      accuracy_peaks = mean(num_segments == 1, na.rm = TRUE),
      N = mean(N),
      eps = mean(eps),
      noise = mean(noise)
      ), by = .(metric, actual_dist)
    ]



table.preds <- table(sim.data.largest$actual_dist, sim.data.largest$dist)
table.preds.prop <- prop.table(table.preds, margin = 1)
addmargins(table.preds, margin = 2)

heatmap.data <- t(table.preds.prop)
heatmap.data <- heatmap.data[c('gamma_flip', 'gamma', 'norm', 'unif'), ]
heatmap.xaxis.lab <- c('Gamma Flip', 'Gamma', 'Normal', 'Uniform')
heatmap.yaxis.lab <- c('Gamma', 'Normal', 'Uniform')

create.heatmap(
  heatmap.data,
  cluster.dimensions = 'none',
  xaxis.lab = heatmap.xaxis.lab,
  yaxis.lab = heatmap.yaxis.lab,
  text.col = 'black',
  cell.text = diag(round(heatmap.data[-1, ], 2)),
  col.pos = c(2,3,4),
  row.pos = c(1,2,3),
  at =  seq(0, 1, length.out = 15),
  colour.scheme = c('white', 'red'),
  colourkey.labels.at = seq(0, 1, length.out = 5),
  colourkey.labels = c('0', '0.25', '0.5', '0.75', '1')
  )

# Boxplot separated by 'Correct' == True
create.boxplot(
  N ~ correct | actual_dist,
  data = sim.data.largest
  )

create.boxplot(
  eps ~ correct | actual_dist,
  data = sim.data.largest
  )

create.boxplot(
  noise ~ correct | actual_dist,
  data = sim.data.largest
  )
