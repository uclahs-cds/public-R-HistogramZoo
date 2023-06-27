library(HistogramZoo);
library(data.table);
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

base.path <- HistogramZoo:::load.config()$root.path;
# sim.file <- file.path(base.path, 'results', 'unimodal_sim_noise', 'processed_results.tsv');
results.folder <- file.path(base.path, 'results', 'unimodal_sim_noise_mle');
plots.folder <- file.path(base.path, 'plots');

sim.tsv.mle.paths <- list.files(results.folder, pattern = 'Unimodal_Sim_MLE.tsv$', full.names = TRUE);
sim.tsv.paths <- list.files(results.folder, pattern = 'Unimodal_Sim.tsv$', full.names = TRUE);

sim.mle.data <- rbindlist(
    lapply(sim.tsv.mle.paths, data.table::fread, sep = '\t'),
    fill = TRUE
    )

sim.data <- rbindlist(
    lapply(sim.tsv.paths, data.table::fread, sep = '\t'),
    fill = TRUE
    )

sim.data <- rbind.data.frame(
  sim.mle.data,
  sim.data
  )

sim.data$seg_length <-  sim.data$end - sim.data$start

write.table(
  sim.data,
  file = file.path(results.folder, 'unimodal_sim_noise.tsv'),
  sep = '\t',
  row.names = FALSE
  )

if (FALSE) {
  sim.data <- read.table(
    file = file.path(results.folder, 'unimodal_sim_noise.tsv'),
    sep = '\t',
    header = TRUE
    )
  }

setDT(sim.data)[, id := .GRP, by = .(N, param, noise, actual_dist, eps, metric)]

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

sim.data.split <- split(sim.data.largest, sim.data.largest$metric)

for (m in names(sim.data.split)) {
  sim.data.metric <- sim.data.split[[m]]

  quantile_p <- c(0, 1/4, 2/4, 3/4, 1)
  sim.data.metric$noise_quartile <- cut(sim.data.metric$noise, quantile(sim.data.metric$noise, probs = quantile_p))
  sim.data.metric$eps_quartile <- cut(sim.data.metric$eps, quantile(sim.data.metric$eps))
  sim.data.metric$N_quartile <- cut(sim.data.metric$N, quantile(sim.data.metric$N))
  levels(sim.data.metric$noise_quartile) <- quantile_p[-1]
  levels(sim.data.metric$eps_quartile) <- quantile_p[-1]
  levels(sim.data.metric$N_quartile) <- quantile_p[-1]

  sim.data.metric$correct <- sim.data.metric$actual_dist == sim.data.metric$dist
  sim.data.metric$correct_fct <- as.factor(sim.data.metric$correct)

  unimodal_accuracy <- sim.data.metric[
      , .(
        accuracy_dist = mean(actual_dist == dist, na.rm = TRUE),
        accuracy_peaks = mean(num_segments == 1, na.rm = TRUE),
        N = mean(N),
        eps = mean(eps),
        noise = mean(noise)
        ), by = .(metric, actual_dist)
      ]

  unimodal_accuracy_filename <- file.path(
    results.folder,
    generate.filename('HZSimulation', paste0(m, '_unimodal_overall_accuracy'), 'tsv')
    )
  cat('Writing unimodal accuracy to: ', unimodal_accuracy_filename, '\n')
  write.table(
    unimodal_accuracy,
    file = unimodal_accuracy_filename,
    sep = '\t',
    row.names = FALSE
    )

  table.preds <- table(sim.data.metric$actual_dist, sim.data.metric$dist)
  table.preds.prop <- prop.table(table.preds, margin = 1)
  # addmargins(table.preds, margin = 2)

  heatmap.data <- t(table.preds.prop)
  heatmap.data <- heatmap.data[c('gamma_flip', 'gamma', 'norm', 'unif'), ]
  heatmap.xaxis.lab <- c('Gamma Flip', 'Gamma', 'Normal', 'Uniform')
  heatmap.yaxis.lab <- c('Gamma', 'Normal', 'Uniform')

  create.heatmap(
    heatmap.data,
    cluster.dimensions = 'none',
    main = m,
    xaxis.lab = heatmap.xaxis.lab,
    yaxis.lab = heatmap.yaxis.lab,
    text.col = 'black',
    cell.text = round(heatmap.data, 2),
    col.pos = rep(rep(1:nrow(heatmap.data)), ncol(heatmap.data)),
    row.pos = unlist(lapply(1:ncol(heatmap.data), function(i) rep(i, nrow(heatmap.data)))),
  #  row.pos = c(1,2,3),
    at =  seq(0, 1, length.out = 15),
    colour.scheme = c('white', 'red'),
    colourkey.labels.at = seq(0, 1, length.out = 5),
    colourkey.labels = c('0', '0.25', '0.5', '0.75', '1'),
    colourkey.cex = 1.5,
    width = 12,
    height = 12,
    resolution = 100,
    filename = file.path(
      plots.folder,
      generate.filename('HZSimulation', paste0(m, '_classification_heatmap'), 'png')
      )
    )

  # Boxplot separated by 'Correct' == True
  n_mle_boxplot <- create.boxplot(
    N ~ correct_fct | actual_dist,
    data = sim.data.metric,
    layout = c(3, 1),
    xlab.label = ''
    )

  eps_mle_boxplot <- create.boxplot(
    eps ~ correct_fct | actual_dist,
    data = sim.data.metric,
    layout = c(3, 1),
    ylab.label = expression(bold('\u03B5')),
    xlab.label = ''
    )

  noise_mle_boxplot <- create.boxplot(
    noise ~ correct_fct | actual_dist,
    data = sim.data.metric,
    layout = c(3, 1),
    xlab.label = ''
    )

  create.multipanelplot(
    list(n_mle_boxplot, eps_mle_boxplot, noise_mle_boxplot),
    main = m,
    layout.width = 3,
    layout.height = 1,
    width = 24,
    height = 8,
    xlab.label = 'Correct Distribution',
    resolution = 100,
    filename = file.path(
      plots.folder,
      generate.filename('HZSimulation', paste0(m, '_mpp_correct_boxplot'), 'png')
      )
  )
}

  # Quartile accuracy of the number of peaks
  quartile.accuracy <- sim.data.metric[
    order(
      actual_dist,
      N_quartile,
      eps_quartile
    ), .(
      accuracy_peaks = mean(num_segments == 1, na.rm = TRUE),
      N = mean(N),
      eps = mean(eps)
      ), by = .(actual_dist, eps_quartile, N_quartile)
    ]

  quartile.accuracy.split <- split(quartile.accuracy, quartile.accuracy$actual_dist)

  quartile.heatmaps <- lapply(names(quartile.accuracy.split), function(d) {
    x <- quartile.accuracy.split[[d]]
    heatmap.data <- dcast(
      x[complete.cases(x)],
      eps_quartile ~ N_quartile, value.var = 'accuracy_peaks'
      )

    heatmap.rownames <- heatmap.data$eps_quartile
    heatmap.data$eps_quartile <- NULL

    heatmap.data <- as.matrix(heatmap.data)

    create.heatmap(
      heatmap.data,
      cluster.dimensions = 'none',
      main = d,
      main.cex = 2,
      text.col = 'black',
      xaxis.lab = c('', colnames(heatmap.data)),
      yaxis.lab = as.character(heatmap.rownames),
      cell.text = round(heatmap.data, 2),
      col.pos = rep(rep(1:nrow(heatmap.data)), ncol(heatmap.data)),
      row.pos = unlist(lapply(1:ncol(heatmap.data), function(i) rep(i, nrow(heatmap.data)))),
      at =  seq(0, 1, length.out = 15),
      colour.scheme = c('white', 'red'),
      colourkey.labels.at = seq(0, 1, length.out = 5),
      colourkey.labels = c('0', '0.25', '0.5', '0.75', '1'),
      colourkey.cex = 1.5
      )
  })

  mpp_heatmap <- file.path(
      plots.folder,
      generate.filename('HZSimulation', 'mpp_correct_segments_heatmap', 'png')
      )
  cat('Saving multipanel heatmap to: ', mpp_heatmap, '\n')
  create.multipanelplot(
    quartile.heatmaps,
    main = 'Number of distribution accuracy',
    layout.width = 3,
    layout.height = 1,
    width = 24,
    height = 8,
    resolution = 100,
    filename = mpp_heatmap
  )

# create.scatterplot(
#   N ~ eps,
#   alpha = 0.1,
#   col = ifelse(sim.data.metric$actual_dist == sim.data.metric$dist, default.colours(2)[1], default.colours(2)[2]),
#   data = sim.data.metric
#   )

# Density lines: Per metric, per distribution, per parameter

# plot(density(sim.data.metric$N[sim.data.metric$correct], from = 0, to = 500))
# lines(density(sim.data.metric$N[! sim.data.metric$correct], from = 0, to = 500), col = 'red')

# Need to plot correct number of segments vs EPS
