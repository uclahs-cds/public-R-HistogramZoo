library(BoutrosLab.plotting.general)
library(HistogramZoo)

df = data.frame(
  "metric_id" = c(0.1, 0.2, 0.3, 0.4),
  "row_id_num_one" = c(1, 1, 1, 1),
  "row_id_num" = c(1, 2, 3, 4),
  "row_id_char" = c("A", "B", "C", "D"),
  "row_id_factor" = factor(LETTERS[1:4], levels = LETTERS[1:4]),
  "colour_id" = c("red", "blue", "red", "green"),
  "start" = c(1,3,5,7),
  "end" = c(2, 4, 6, 8),
  "start_id" = c(1,3,5,7),
  "end_id" = c(2, 4, 6, 8)
)

create_trackplot(
  track_data = df,
  row_id = "row_id_num",
  metric_id = "metric_id",
  colour_scheme = c("white", "black")
)

create_trackplot(
  track_data = df,
  row_id = "row_id_char",
  metric_id = "metric_id",
  colour_scheme = c("white", "black")
)

create_trackplot(
  track_data = df,
  row_id = "row_id_factor",
  metric_id = "metric_id",
  colour_scheme = c("white", "black")
)

create_trackplot(
  track_data = df,
  row_id = "row_id_factor",
  metric_id = "metric_id",
  colour_id = "colour_id",
  colour_scheme = c("white", "black")
)

create_trackplot(
  track_data = df,
  row_id = "row_id_num_one",
  metric_id = "metric_id",
  colour_id = "colour_id"
)

create_trackplot(
  track_data = df,
  row_id = "row_id_num_one",
  metric_id = "metric_id",
  colour_id = "colour_id",
  start_id = "start_id",
  end_id = "end_id"
)
