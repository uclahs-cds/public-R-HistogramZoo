library(ConsensusPeaks)
# Instead of using MLE, use the metrics we are interested in:
# Jaccard or histogram intersection

metric = "jaccard"

hist.optim.plot <- function(x, metric = "jaccard", dist = c('norm', 'gamma')) {
   if(dist == "gamma" && min(x) < 0) {
      x = x - min(x)
   }
   bin.data = obs.to.int.hist(x, as.data.frame = TRUE, add.zero.endpoints = FALSE)

   dist.optim = fit.distributions.optim(bin.data, metric = metric, truncated = truncated)

   dens <- dnorm(bin.data$x, mean = dist.optim$norm$par[1], sd = dist.optim$norm$par[2]) * sum(bin.data$Freq)

   plot(bin.data$x, dens,
        col = scales::alpha("pink", 0.7), type = "h", lwd = 5, lend = 1, ylim = c(0, max(dens, bin.data$Freq)),
        main = sprintf("%s optimization for normal distribution", metric))
   lines(bin.data$x, bin.data$Freq, col = scales::alpha("lightblue", 0.9), type = "h", lwd = 5, lend = 1)
   lines(bin.data$x, pmin(a = bin.data$Freq, b = dens, na.rm = T),
         col = scales::alpha("lightgreen", 0.9), type = "h", lwd = 5, lend = 1)
   legend("topleft",
          legend = c("Observed", "Model", "Intersect"),
          fill = c(scales::alpha("lightblue", 0.9), scales::alpha("pink", 0.7), scales::alpha("lightgreen", 0.9)),
          cex = 0.5,
          inset = .05)
   mtext(sprintf("%s index: %.3f", metric, round(1 - dist.optim$norm$value, 3)))
}

set.seed(13)
mu = 0
sigma = 50
x <- rnorm(1e4, mean = mu, sd = sigma)
bin.data = obs.to.int.hist(x - min(x), as.data.frame = TRUE, add.zero.endpoints = FALSE)

dist.optim = fit.distributions.optim(bin.data, metric = metric, truncated = FALSE)
# intersect.optim = fit.distributions.optim(bin.data, metric = "intersection")

dens.norm <- dnorm(bin.data$x, mean = dist.optim$norm$par[1], sd = dist.optim$norm$par[2]) * sum(bin.data$Freq)

plot(bin.data$x, dens.norm,
     col = scales::alpha("pink", 0.7), type = "h", lwd = 5, lend = 1, ylim = c(0, max(dens.norm, bin.data$Freq)),
     main = sprintf("%s optimization for normal distribution", metric))
lines(bin.data$x, bin.data$Freq, col = scales::alpha("lightblue", 0.9), type = "h", lwd = 5, lend = 1)
lines(bin.data$x, pmin(a = bin.data$Freq, b = dens.norm, na.rm = T),
      col = scales::alpha("lightgreen", 0.9), type = "h", lwd = 5, lend = 1)
legend("topleft",
       legend = c("Observed", "Model", "Intersect"),
       fill = c(scales::alpha("lightblue", 0.9), scales::alpha("pink", 0.7), scales::alpha("lightgreen", 0.9)),
       cex = 0.5,
       inset = .05)
mtext(sprintf("%s index: %.3f", metric, round(1 - dist.optim$norm$value, 3)))

dens.gamma <- dgamma(bin.data$x, shape = dist.optim$gamma$par[1], rate = dist.optim$gamma$par[2]) * sum(bin.data$Freq)

plot(bin.data$x, dens.gamma,
     col = scales::alpha("pink", 0.7), type = "h", lwd = 5, lend = 1, ylim = c(0, max(dens.gamma, bin.data$Freq)),
     main = sprintf("%s optimization for gamma distribution", metric))
lines(bin.data$x, bin.data$Freq, col = scales::alpha("lightblue", 0.9), type = "h", lwd = 5, lend = 1)
lines(bin.data$x, pmin(a = bin.data$Freq, b = dens.gamma, na.rm = T),
      col = scales::alpha("lightgreen", 0.9), type = "h", lwd = 5, lend = 1)
legend("topleft",
       legend = c("Observed", "Model", "Intersect"),
       fill = c(scales::alpha("lightblue", 0.9), scales::alpha("pink", 0.7), scales::alpha("lightgreen", 0.9)),
       cex = 0.5,
       inset = .05)
mtext(sprintf("%s index: %.3f", metric, round(1 - dist.optim$gamma$value, 3)))

# Use a unif distribution
set.seed(13)
x <- runif(1e4, min = 0, max = 100)
bin.data = obs.to.int.hist(100 * (x - min(x)) / (max(x) - min(x)), as.data.frame = TRUE, add.zero.endpoints = FALSE)

dist.optim = fit.distributions.optim(bin.data, metric = metric, truncated = TRUE)

dens.norm <- dtnorm(bin.data$x, mean = dist.optim$norm$par[1], sd = dist.optim$norm$par[2],
                   a = min(bin.data$x) - 1e-10, b = max(bin.data$x) + 1e-10) * sum(bin.data$Freq)

plot(bin.data$x, dens.norm,
     col = scales::alpha("pink", 0.7), type = "h", lwd = 5, lend = 1, ylim = c(0, max(dens.norm, bin.data$Freq)),
     main = sprintf("%s optimization for normal distribution", metric))
lines(bin.data$x, bin.data$Freq, col = scales::alpha("lightblue", 0.9), type = "h", lwd = 5, lend = 1)
lines(bin.data$x, pmin(a = bin.data$Freq, b = dens.norm, na.rm = T),
      col = scales::alpha("lightgreen", 0.9), type = "h", lwd = 5, lend = 1)
legend("topleft",
       legend = c("Observed", "Model", "Intersect"),
       fill = c(scales::alpha("lightblue", 0.9), scales::alpha("pink", 0.7), scales::alpha("lightgreen", 0.9)),
       cex = 0.5,
       inset = .05)
mtext(sprintf("%s index: %.3f", metric, round(1 - dist.optim$norm$value, 3)))

dens.gamma <- dtgamma(bin.data$x, shape = dist.optim$gamma$par[1], rate = dist.optim$gamma$par[2],
                      a = min(bin.data$x) - 1e-10, b = max(bin.data$x) + 1e-10) * sum(bin.data$Freq)

plot(bin.data$x, dens.gamma,
     col = scales::alpha("pink", 0.7), type = "h", lwd = 5, lend = 1, ylim = c(0, max(dens.gamma, bin.data$Freq)),
     main = sprintf("%s optimization for gamma distribution", metric))
lines(bin.data$x, bin.data$Freq, col = scales::alpha("lightblue", 0.9), type = "h", lwd = 5, lend = 1)
lines(bin.data$x, pmin(a = bin.data$Freq, b = dens.gamma, na.rm = T),
      col = scales::alpha("lightgreen", 0.9), type = "h", lwd = 5, lend = 1)
legend("topleft",
       legend = c("Observed", "Model", "Intersect"),
       fill = c(scales::alpha("lightblue", 0.9), scales::alpha("pink", 0.7), scales::alpha("lightgreen", 0.9)),
       cex = 0.5,
       inset = .05)
mtext(sprintf("%s index: %.3f", metric, round(1 - dist.optim$gamma$value, 3)))

# Use a unif distribution, without truncation
set.seed(13)
x <- runif(1e4, min = 0, max = 100)
bin.data = obs.to.int.hist(100 * (x - min(x)) / (max(x) - min(x)), as.data.frame = TRUE, add.zero.endpoints = FALSE)

dist.optim = fit.distributions.optim(bin.data, metric = metric, truncated = FALSE)

dens.norm <- dnorm(bin.data$x, mean = dist.optim$norm$par[1], sd = dist.optim$norm$par[2]) * sum(bin.data$Freq)

dens.gamma <- dtgamma(bin.data$x, shape = dist.optim$gamma$par[1], rate = dist.optim$gamma$par[2]) * sum(bin.data$Freq)

plot(bin.data$x, dens.norm,
     col = scales::alpha("pink", 0.7), type = "h", lwd = 5, lend = 1, ylim = c(0, max(dens.norm, bin.data$Freq)),
     main = sprintf("%s optimization for normal distribution", metric))
lines(bin.data$x, bin.data$Freq, col = scales::alpha("lightblue", 0.9), type = "h", lwd = 5, lend = 1)
lines(bin.data$x, pmin(a = bin.data$Freq, b = dens.norm, na.rm = T),
      col = scales::alpha("lightgreen", 0.9), type = "h", lwd = 5, lend = 1)
legend("topleft",
       legend = c("Observed", "Model", "Intersect"),
       fill = c(scales::alpha("lightblue", 0.9), scales::alpha("pink", 0.7), scales::alpha("lightgreen", 0.9)),
       cex = 0.5,
       inset = .05)
mtext(sprintf("%s index: %.3f", metric, round(1 - dist.optim$norm$value, 3)))

plot(bin.data$x, dens.gamma,
     col = scales::alpha("pink", 0.7), type = "h", lwd = 5, lend = 1, ylim = c(0, max(dens.gamma, bin.data$Freq)),
     main = sprintf("%s optimization for gamma distribution", metric))
lines(bin.data$x, bin.data$Freq, col = scales::alpha("lightblue", 0.9), type = "h", lwd = 5, lend = 1)
lines(bin.data$x, pmin(a = bin.data$Freq, b = dens.gamma, na.rm = T),
      col = scales::alpha("lightgreen", 0.9), type = "h", lwd = 5, lend = 1)
legend("topleft",
       legend = c("Observed", "Model", "Intersect"),
       fill = c(scales::alpha("lightblue", 0.9), scales::alpha("pink", 0.7), scales::alpha("lightgreen", 0.9)),
       cex = 0.5,
       inset = .05)
mtext(sprintf("%s index: %.3f", metric, round(1 - dist.optim$gamma$value, 3)))
