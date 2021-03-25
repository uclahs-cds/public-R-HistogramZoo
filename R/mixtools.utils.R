dnorm_mixture <- function(x, mixfit) {
  Reduce('+', lapply(1:length(mixfit$lambda), function(i) mixfit$lambda[i] * dnorm(x, mixfit$mu[i], mixfit$sigma[i])))
}
