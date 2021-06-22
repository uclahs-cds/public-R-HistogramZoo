library(dplyr)
library(tidyr)
library(BoutrosLab.plotting.general)
devtools::load_all(".")

data <- readRDS("rds/foxa1.rds")
attach(data)

# Only look at second peak
filenames <- NULL
metrics = c("intersect", "jaccard")
I = c(2, 4)
for(metric in metrics) {
  for(i in I) {
    seg.start = GenomicRanges::start(seg.gr)[i]
    seg.end = GenomicRanges::end(seg.gr)[i]
    #
    x = peak.counts[peak.counts >= seg.start & peak.counts <= seg.end]
    x.adj = x - min(x) + 1
    bin.data = obs.to.int.hist(x.adj, as.data.frame = TRUE, add.zero.endpoints = FALSE)
    dist.optim = fit.distributions.optim(bin.data, metric = metric, truncated = FALSE)

    filename <- sprintf("plots/foxa1_%s_seg_%s.png", metric, i)
    filenames <- c(filenames, filename)
    png(filename = filename, width = 8, height = 8, units = "in", res = 500)
    par(oma = c(0, 0, 2, 0))
    hist.optim.plot.all(bin.data, dist.optim, metric = metric)
    mtext(sprintf("FOXA1 Segment %d: %s", i, metric), side = 3, line = 0, outer = TRUE, cex = 2)
    dev.off()

    # Fit it on the largest uniform segment
    max.unif.results <- find.uniform.segment(bin.data$Freq, threshold = .75, step.size = 5)
    unif.segment = unlist(max.unif.results[c('a', 'b')])
    unif.segment.adj =  unif.segment + min(bin.data$x)
    x.subset <- x[x >= unif.segment.adj[1] & x <= unif.segment.adj[2]]
    bin.data.subset <- bin.data[bin.data$x >= unif.segment.adj[1] & bin.data$x <= unif.segment.adj[2],]
    bin.data.subset$x <- bin.data.subset$x - unif.segment.adj[1] + 1
    dist.optim.subset = fit.distributions.optim(bin.data.subset, metric = metric, truncated = FALSE)

    filename <- sprintf("plots/foxa1_%s_uniform_seg_%s.png", metric, i)
    filenames <- c(filenames, filename)
    png(filename = filename, width = 8, height = 8, units = "in", res = 500)
    par(oma = c(0, 0, 2, 0))
    hist.optim.plot.all(bin.data, dist.optim.subset, metric = metric, bounds = unif.segment.adj)
    mtext(sprintf("FOXA1 Max Uniform Segment %d: %s", i, metric), side = 3, line = 0, outer = TRUE, cex = 2)
    dev.off()
  }
}

cat("convert", paste0(filenames, collapse = " "), "plots/foxa1_histogram_optim.pdf")
