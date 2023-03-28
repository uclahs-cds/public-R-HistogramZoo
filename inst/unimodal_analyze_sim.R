library(HistogramZoo);
library(data.table);
library(BoutrosLab.plotting.general);

base.path <- HistogramZoo:::load.config()$root.path;
sim.file <- file.path(base.path, 'results', 'unimodal_sim_noise', 'processed_results.tsv');
sim.data <- fread(sim.file, sep = '\t')

sim.data$seg_length <-  sim.data$end - sim.data$start

# Look at segments > 1
# sim.data[, .N, by = .(id, metric)][N > 1, .N, by = .(id)]$id

# Only keep the largest segment
sim.data.largest <- sim.data[sim.data[, .I[which.max(seg_length)], by=.(id, metric)]$V1]

# Overall Accuracy
sim.data.largest[, .(correct = actual_dist == dist), by = .(id, metric, N, eps)
    ][
    , .(
      accuracy = mean(correct, na.rm = TRUE),
      N = mean(N),
      eps = mean(eps)
      ), by = .(metric)
    ]

# Segmentation level (ignore metric)
id.group <- sim.data.largest[, .(
  correct_segment = max(num_segments) == 1,
  time = mean(timing),
  N = floor(N),
  eps = max(eps)
  ), by = .(id, actual_dist)]

# TODO: Look more at length of other segments
id.group[, .(
  accuracy_segment = mean(correct_segment)
  ), by = actual_dist]

# summary(glm(correct_segment ~ N + eps + actual_dist, data = id.group))

# Metric level grouping
metric.group <- sim.data.largest[, .(
  N = floor(N),
  eps = max(eps),
  correct_dist = actual_dist == dist
  ), by = .(id, metric, actual_dist, dist, dist_param1, dist_param2, param)]

metric.group[, .(
  acc = mean(correct_dist)
  ), by = .(metric, actual_dist)][order(-acc), ]


# If correct, how far away are parameter estimates?
# Parameters are on histogram bin scale
metric.group[
  actual_dist == dist & actual_dist == 'norm',
  .(mean = mean(dist_param1),
    sd = mean(dist_param2)),
  by = .(metric)]

summarize_results(
  segment_and_fit(observations_to_histogram(rnorm(1e4, mean = -10))))
