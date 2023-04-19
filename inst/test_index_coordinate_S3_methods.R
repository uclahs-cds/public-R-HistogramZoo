

library(HistogramZoo)
library(BoutrosLab.plotting.general)

# Preamble ----------------------------------------------------------------
# Writing a set of tests to test S3 methods for fitting histograms
# Updated functions

# Coordinate based functions
# 1. fit_distributions
# 2. fit_uniform
# 3. identify_uniform_segment

# Index based functions
# 4. find_local_optima
# 5. FTC
# 6. meaningful_gaps_local, find_all_meaningful_gap

# Plotting function -------------------------------------------------------

distribution_colours <- c(
  "norm" = "darkorange",
  "gamma" = "chartreuse4",
  "gamma_flip" = "chartreuse3",
  "unif" = "darkorchid4"
)

metric_lty <- c(
  "mle" = 1,
  "jaccard" = 2,
  "intersection" = 3,
  "chisq" = 4,
  "ks" = 5,
  "mse" = 6
)

plot_fit_distributions <- function(res, x, empirical_data){
  
  plotting_data <- do.call(rbind.data.frame, lapply(res, function(res_i){
    data.frame(
      "x" = x,
      "y" = res_i$dens(x, mpar = res_i$par, scale = T),
      "metric" = res_i[['metric']],
      "distribution" = res_i[['dist']],
      "tag" = paste0(res_i[['metric']], "-", res_i[['dist']])
    )
  }))
  
  plt <- create.scatterplot(
    y ~ x,
    plotting_data,
    # plotting data
    groups = plotting_data$tag,
    type = "l",
    col = distribution_colours[unique(plotting_data$distribution)],
    lty = metric_lty[unique(plotting_data$metric)],
    # Axes
    xlab.cex = 1,
    xlab.label = "x",
    ylab.cex = 1,
    ylab.label = "Data",
    xaxis.cex = 1,
    yaxis.cex = 1,
    xlimits = c(min(x), max(x)),
    # Adding points
    add.points = T,
    points.x = x,
    points.y = empirical_data,
    points.col = "red"
  )
  
  return(plt)
}

find_midpoint <- function(start, end){
  (start + end)/2
}

# Gamma testing -----------------------------------------------------------

# Gamma histogram
gamma_data <- rgamma(100000, shape = 20, rate = 2, shift = 100)
gamma_density <- dgamma(1:200, shape = 20, rate = 2, shift = 100)
gamma_histogram <- observations_to_histogram(gamma_data)
res <- fit_distributions(gamma_histogram, dist = "gamma", truncated = T)
plot_fit_distributions(
  res = res, 
  x = 1:200, 
  gamma_density
)

# Gamma test
gamma_data <- dgamma(1:200, shape = 20, rate = 2, shift = 100)
tgamma_data <- dtgamma(1:200, shape = 20, shift = 100, rate = 2, a = 10, b = 200)
res <- HistogramZoo:::fit_distributions_helper(
  gamma_data[100:200], 
  dist = "gamma", 
  metric = "jaccard", 
  truncated = T,
  interval_start = 100:200 - 0.5,
  interval_end = 100:200 + 0.5,
  interval_midpoint = 100:200
)

plot_fit_distributions(
  res = res, 
  x = 100:200, 
  gamma_data[100:200]
)

plot(x = 1:200, y = gamma_data, col = "red")
points(x = 1:200, y = tgamma_data, col = "blue")

# Gamma flip example
gamma_flip_data <- rgamma_flip(100000, shape = 20, rate = 2, offset = 200)
gamma_flip_density <- dgamma_flip(1:200, shape = 20, rate = 2, offset = 200)
gamma_flip_histogram <- observations_to_histogram(gamma_flip_data)
res <- fit_distributions(gamma_flip_histogram, dist = "gamma_flip", truncated = F)
plot_fit_distributions(
  res = res, 
  x = 1:200,
  gamma_flip_density
)

# Gamma flip test
gamma_flip_data <- dgamma_flip(1:200, shape = 20, rate = 2, offset = 200)
tgamma_flip_data <- dtgamma_flip(1:200, shape = 20, rate = 2, offset = 200, a = 0, b = Inf)
res <- HistogramZoo:::fit_distributions_helper(
  gamma_flip_data[100:200], 
  dist = "gamma_flip", 
  truncated = T,
  interval_start = 100:200 - 0.5,
  interval_end = 100:200 + 0.5,
  interval_midpoint = 100:200
)

plot_fit_distributions(
  res = res, 
  x = 1:200,
  gamma_data
)

plot(x = 1:200, y = gamma_flip_data, col = "red")
points(x = 1:200, y = tgamma_flip_data, col = "blue")

