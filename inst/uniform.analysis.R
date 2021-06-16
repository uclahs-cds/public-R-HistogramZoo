library(dplyr)
library(tidyr)
library(BoutrosLab.plotting.general)

data <- readRDS("rds/foxa1.rds")
attach(data)

hist(peak.counts, breaks = seq(1, max(peak.counts)))

seg.start = GenomicRanges::start(seg.gr)[2]
seg.end = GenomicRanges::end(seg.gr)[2]

x = peak.counts[peak.counts >= seg.start & peak.counts <= seg.end]


x.adjusted <- (x - seg.start) + 1e-10
x.range = seg.start:seg.end
x.range.adjusted <- (x - seg.start) + 1e-10

hist(x.adjusted, breaks = seq(0, floor(max(x.adjusted))))

# Refit models
mod = fit.continuous.distributions(
  x = x.adjusted,
  seg.start = seg.start,
  seg.end = seg.end,
  fit.mixtures = c('unif', 'tnorm', 'tgamma', 'tgamma_flip'),
  max.iterations = 500)

fits = extract.distribution.parameters(
  mod = mod,
  x = x.adjusted)


res <- lapply(seq(105, 155, by = 10), function(left.adj) {
  print(left.adj)
  forced.unif.x <- x.adjusted[x.adjusted >= left.adj]
  mod.unif <- fit.continuous.distributions(
    x = forced.unif.x,
    seg.start = seg.start + left.adj,
    seg.end = seg.end,
    fit.mixtures = c('unif', 'tnorm', 'tgamma'),
    max.iterations = 500)

  fits.unif = extract.distribution.parameters(
    mod = mod.unif,
    x = forced.unif.x)
  fits.unif$left.adj <- left.adj
  fits.unif
})


res.merge <- do.call(rbind.data.frame, res)

res.merge %>%
  group_by(left.adj) %>%
  slice_max(jc, n = 1)


# Uniform segment idenification
identify.uniform.segments(seg.gr, peak.counts, trim.peak.threshold = 0.5, short.peak.threshold = 10)

res <- find.uniform.segment(x, step.size = 2)
res$index <- 1:nrow(res)
create.segplot(formula = index ~ a + b, data = res)

plot.data = data.frame(peak.counts = seq(1, max(peak.counts)))
plot.data = merge(plot.data, as.data.frame(table(peak.counts)), all.x = TRUE)
plot.data$Freq[is.na(plot.data$Freq)] <- 0

plot.data$unif.jaccard = 0

segments <- list(c(749, 1078), c(1317, 3440), c(3341, 3740), c(3741, 4414))
unif.results = do.call(rbind.data.frame, lapply(seq_along(segments), function(i) {
  seg.start = segments[[i]][1]
  seg.end = segments[[i]][2]

  x = plot.data[plot.data$peak.counts >= seg.start & plot.data$peak.counts <= seg.end, "Freq"]
  x.adjusted <- (x - seg.start)
#  x.range = seg.start:seg.end
#  x.range.adjusted <- (x - seg.start + 1)

  # x.table <- table(x.adjusted)

  unif.results <- find.uniform.segment(x.adjusted, step.size = 20, threshold = .85)
  unif.results[, c('a', 'b')] <- unif.results[, c('a', 'b')] + seg.start - 1
  unif.results$n <- unif.results$b - unif.results$a + 1
  unif.results$L <- length(x.adjusted)

  unif.results %>%
    mutate(neg.jaccard = -jaccard.index, seg = i) %>%
    pivot_longer(c(entropy, mse, neg.jaccard)) %>%
    group_by(name) %>%
    slice_min(value, with_ties = FALSE)
}))

#jaccard.heatmap =

unif.jaccard = unlist(apply(unif.results[unif.results$name == "neg.jaccard", ], 1, function(r) {
  seq(r[["a"]], r[["b"]])
}))

plot.data$unif.jaccard[unif.jaccard] <- 1

# rownames(plot.data)  = plot.data$peak.counts
main = create.scatterplot(
  Freq ~ peak.counts,
  plot.data,
#   xaxis.cex = 0,
  xlab.cex = 0,
  ylab.cex = 1,
#  xaxis.tck = 0,
  yaxis.cex = 1,
  ylab.label = "Coverage (at BP resolution)",
  main.cex = 0,
  type = "a",
  add.points = T,
  points.x = p,
  points.y = plot.data$Freq[p],
  points.pch = 19,
  points.col = 'red'
)

jc.heatmap = create.heatmap(
  t(as.matrix(plot.data$unif.jaccard)),
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

res2 <- find.uniform.segment(table(peak.counts), step.size = 10, threshold = 0.6)
res2$index <- 1:nrow(res2)
create.segplot(formula = index ~ a + b, data = res2)

res2 %>%
  mutate(neg.jaccard = -jaccard.index) %>%
  pivot_longer(c(mse, neg.jaccard)) %>%
  group_by(name) %>%
  slice_min(value, with_ties = FALSE)
