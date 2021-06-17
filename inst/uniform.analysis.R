library(dplyr)
library(tidyr)
library(BoutrosLab.plotting.general)

data <- readRDS("rds/foxa1.rds")
devtools::load_all(".")
attach(data)

hist(peak.counts, breaks = seq(1, max(peak.counts)))

seg.start = GenomicRanges::start(seg.gr)[2]
seg.end = GenomicRanges::end(seg.gr)[2]

x = peak.counts[peak.counts >= seg.start & peak.counts <= seg.end]


#x.adjusted <- (x - seg.start) + 1e-10
#x.range = seg.start:seg.end
#x.range.adjusted <- (x - seg.start) + 1e-10

x.df <- as.data.frame(table(x))
x.df$x <- as.numeric(as.character(x.df$x))
unif.data = x.df
# unif.data = data.frame(peak.counts = seq(1, max(x)))
# unif.data = merge(unif.data, x.df, all.x = TRUE)

# hist(x.adjusted, breaks = seq(0, floor(max(x.adjusted))))
res <- find.uniform.segment(x.df$Freq, threshold = .75, step.size = 5)
unif.segment = unlist(res[c('a', 'b')])  + min(x.df$x)
points.y = x.df$Freq[unlist(res[c('a', 'b')])]
main = create.scatterplot(
  Freq ~ x,
  x.df,
  #   xaxis.cex = 0,
  xlab.cex = 0,
  ylab.cex = 1,
  #  xaxis.tck = 0,
  yaxis.cex = 1,
  xlimits = c(min(x.df$x), max(x.df$x)),
  ylab.label = "Coverage (at BP resolution)",
  main.cex = 0,
  type = "a",
  add.points = T,
  points.x = unif.segment,
  points.y = points.y,
  points.pch = 19,
  points.col = 'red'
)

unif.data$unif <- unif.data$x >= unif.segment[1] & unif.data$x <= unif.segment[2]
jc.heatmap = create.heatmap(
  t(as.matrix(unif.data$unif)),
  clustering.method = 'none',
  yaxis.tck = 0,
  total.colours = 3,
  colour.scheme = c("white", "black"),
  print.colour.key = F
)

create.multipanelplot(
  plot.objects = list(main, jc.heatmap),
  plot.objects.heights = c(7, 0.5),
  y.spacing = -1,
  ylab.label = "Start",
  xlab.label = "End",
  ylab.cex = 2,
  xlab.cex = 2,
  height = 8,
  width = 12,
  filename = "plots/test_jc.png"
)
