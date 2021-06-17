library(ConsensusPeaks)
# Instead of using MLE, use the metrics we are interested in
# Jaccard or histogram intersection
# Let's see how it works for fitting a normal distribution

set.seed(13)
mu = 0
sigma = 50
x <- rnorm(1e4, mean = mu, sd = sigma)
bin.data = obs.to.int.hist(x, as.data.frame = TRUE, add.zero.endpoints = FALSE)

dens <- dnorm(bin.data$x, mean = mu, sd = sigma) * sum(bin.data$Freq)
plot(bin.data$x, dens, col = scales::alpha("pink", 0.7), type = "h", lwd = 5, lend = 1, ylim = c(0, max(dens, bin.data$Freq)))
lines(bin.data$x, bin.data$Freq, col = scales::alpha("lightblue", 0.9), type = "h", lwd = 5, lend = 1)
lines(bin.data$x, pmin(a = bin.data$Freq, b = dens, na.rm = T),
      col = scales::alpha("lightgreen", 0.9), type = "h", lwd = 5, lend = 1)
legend("topleft",
       legend = c("Observed", "Model", "Intersect"),
       fill = c(scales::alpha("lightblue", 0.9), scales::alpha("pink", 0.7), scales::alpha("lightgreen", 0.9)),
       cex = 0.5,
       inset = .05)


# Histogram comparison
# https://stats.stackexchange.com/a/151362/97417
normal.hist.optim <- function(data, params, metric = c("jaccard", "intersect")) {
  metric = match.arg(metric)
  dens <- dnorm(data$x, mean = params[1], sd = params[2]) * sum(data$Freq)
  overlap = pmin(a = data$Freq, b = dens, na.rm = T)
  if(metric == "jaccard") {
    union = pmax(a = data$Freq, b = dens, na.rm = T)
    rtn = sum(overlap)/sum(union)
  } else if (metric == "intersect") {
    rtn = sum(overlap) / sum(data$Freq)
  }

  rtn
}

res <- optim(par = c(0, 1),
      fn = normal.hist.optim,
      method = "L-BFGS-B",
      data = bin.data,
      control = list(fnscale = -1),
      lower = c(-Inf, 0),
      metric = "intersect"
      )

# Use a unif distribution
set.seed(13)
x <- runif(1e4, min = 0, max = 100)
bin.data = obs.to.int.hist(x, as.data.frame = TRUE, add.zero.endpoints = FALSE)
res <- optim(par = c(0, 1),
             fn = normal.hist.optim,
             method = "L-BFGS-B",
             data = bin.data,
             control = list(fnscale = -1),
             lower = c(-Inf, 0),
             metric = "intersect")
dens <- dnorm(bin.data$x, mean = res$par[1], sd = res$par[2]) * sum(bin.data$Freq)
plot(bin.data$x, dens, col = scales::alpha("pink", 0.7), type = "h", lwd = 5, lend = 1, ylim = c(0, max(dens, bin.data$Freq)))
lines(bin.data$x, bin.data$Freq, col = scales::alpha("lightblue", 0.9), type = "h", lwd = 5, lend = 1)
lines(bin.data$x, pmin(a = bin.data$Freq, b = dens, na.rm = T),
      col = scales::alpha("lightgreen", 0.9), type = "h", lwd = 5, lend = 1)
legend("bottomleft",
       legend = c("Observed", "Model", "Intersect"),
       fill = c(scales::alpha("lightblue", 0.9), scales::alpha("pink", 0.7), scales::alpha("lightgreen", 0.9)),
       cex = 0.5,
       inset = 0.05
       )
