library(dplyr)
library(tidyr)
library(BoutrosLab.plotting.general)
devtools::load_all(".")

data <- readRDS("rds/foxa1.rds")
attach(data)

# Only look at second peak
seg.start = GenomicRanges::start(seg.gr)[2]
seg.end = GenomicRanges::end(seg.gr)[2]

x = peak.counts[peak.counts >= seg.start & peak.counts <= seg.end]

# Table the data
x.df <- as.data.frame(table(x))
x.df$x <- as.numeric(as.character(x.df$x))
unif.data = x.df

res <- find.uniform.segment(x.df$Freq, threshold = .75, step.size = 5)
unif.segment = unlist(res[c('a', 'b')])
unif.segment.adj =  unif.segment + min(x.df$x)
points.y = x.df$Freq[unif.segment]
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
  points.x = unif.segment.adj,
  points.y = points.y,
  points.pch = 19,
  points.col = 'red'
)

unif.data$unif <- unif.data$x >= unif.segment.adj[1] & unif.data$x <= unif.segment.adj[2]
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

# Fit the model on the original data
x.adjusted <- (x - unif.segment[1]) + 1e-10
mod = fit.continuous.distributions(
  x = x.adjusted,
  seg.start = seg.start,
  seg.end = seg.end,
  fit.mixtures = c('unif', 'tnorm', 'tgamma', 'tgamma_flip'),
  max.iterations = 500)

fits = extract.distribution.parameters(
  mod = mod,
  x = x.adjusted)

# Fit the models on the unif segment
x.segment <- x[x >= unif.segment.adj[1] & x <= unif.segment.adj[2]]
x.segment.adjusted <- (x.segment - unif.segment[1]) + 1e-10
mod.unif = fit.continuous.distributions(
  x = x.segment.adjusted,
  seg.start =  unif.segment[1],
  seg.end = unif.segment[2],
  fit.mixtures = c('unif', 'tnorm', 'tgamma', 'tgamma_flip'),
  max.iterations = 500)

fits.unif = extract.distribution.parameters(
  mod = mod.unif,
  x = x.segment.adjusted)
