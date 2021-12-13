hist.optim.plot <- function(bin.data, dist.fit = NULL, metric = "jaccard", dist = c('norm', 'gamma', 'unif'), truncated = FALSE, bounds = NULL, legend = TRUE) {
  if(!is.null(bounds)) {
    bin.data.subset <- bin.data[bin.data$x >= bounds[1] & bin.data$x <= bounds[2],]
    bin.data.subset$x <- bin.data.subset$x - bounds[1] + 1
  } else {
    bin.data.subset <- bin.data
  }
  if(is.null(dist.fit)) {
    dist.fit = fit.distributions.optim(bin.data.subset, metric = metric, truncated = FALSE)
  }
  dens = dist.fit[[dist]]$dens()

  if(!is.null(bounds)) {
    bin.data.subset$x <- bin.data.subset$x + bounds[1] + 1
  }

  plot(bin.data.subset$x, dens,
       col = scales::alpha("pink", 0.7), type = "h", lwd = 5, lend = 1, xlim = c(0, max(bin.data$x)), ylim = c(0, max(bin.data$Freq) + round(sd(bin.data$Freq))),
       main = sprintf("%s optimization for %s distribution", metric, dist),
       xlab = "x", ylab = "dens")
  lines(bin.data$x, bin.data$Freq, col = scales::alpha("lightblue", 0.9), type = "h", lwd = 5, lend = 1)
  lines(bin.data.subset$x, pmin(a = bin.data.subset$Freq, b = dens, na.rm = T),
        col = scales::alpha("lightgreen", 0.9), type = "h", lwd = 5, lend = 1)
  if(legend) {
    legend("topleft",
           legend = c("Observed", "Model", "Intersect"),
           fill = c(scales::alpha("lightblue", 0.9), scales::alpha("pink", 0.7), scales::alpha("lightgreen", 0.9)),
           cex = 0.5,
           inset = .05)
  }

  if(metric == "ks") {
    metric.value <- dist.fit[[dist]]$value
  } else {
    metric.value <- 1 - dist.fit[[dist]]$value
  }

  mtext(sprintf("%s index: %.3f", metric, metric.value))
}

hist.optim.plot.all <- function(bin.data, dist.fit = NULL, metric = "jaccard", truncated = FALSE, bounds = NULL) {
  par(mfrow = c(2,2))
  hist.optim.plot(bin.data, dist.fit, dist = "norm", metric = metric, bounds = bounds, legend = FALSE)
  hist.optim.plot(bin.data, dist.fit, dist = "gamma", metric = metric, bounds = bounds, legend = FALSE)
  hist.optim.plot(bin.data, dist.fit, dist = "unif", metric = metric, bounds = bounds, legend = FALSE)
  plot(1, type = "n", xlab = "", ylab = "", axes = FALSE)
  legend("topleft",
         legend = c("Observed", "Model", "Intersect"),
         fill = c(scales::alpha("lightblue", 0.9), scales::alpha("pink", 0.7), scales::alpha("lightgreen", 0.9)),
         cex = 0.75,
         inset = .05,
         horiz = TRUE)
}
