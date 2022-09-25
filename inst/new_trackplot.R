
library(BoutrosLab.plotting.general)

df = data.frame(
  "start" = c(1:5),
  "end" = c(6:10),
  "row" = c(1:5),
  "metric" = c(1:5)
)

create.scatterplot(
  formula = start ~ end,
  data = df,
  col = "transparent",
  add.rectangle = T,
  xleft.rectangle = df$start,
  ybottom.rectangle = df$row - 0.5,
  xright.rectangle = df$end,
  ytop.rectangle = df$row + 0.5,
  col.rectangle = c("pink", "red"),
  alpha.rectangle = 1,
)