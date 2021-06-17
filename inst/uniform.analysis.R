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


x.adjusted <- (x - seg.start) + 1e-10
x.range = seg.start:seg.end
x.range.adjusted <- (x - seg.start) + 1e-10

x.df <- as.data.frame(table(x))
x.df$x <- as.numeric(as.character(x.df$x))
# hist(x.adjusted, breaks = seq(0, floor(max(x.adjusted))))
res <- find.uniform.segment(x.df$Freq, threshold = .75, step.size = 5)
points.x = unlist(res[c('a', 'b')]) + min(x.df$x)
points.y = x.df$Freq[unlist(res[c('a', 'b')])]
create.scatterplot(
  Freq ~ x,
  x.df,
  #   xaxis.cex = 0,
  xlab.cex = 0,
  ylab.cex = 1,
  #  xaxis.tck = 0,
  yaxis.cex = 1,
  ylab.label = "Coverage (at BP resolution)",
  main.cex = 0,
  type = "a",
  add.points = T,
  points.x = points.x,
  points.y = points.y,
  points.pch = 19,
  points.col = 'red'
)
