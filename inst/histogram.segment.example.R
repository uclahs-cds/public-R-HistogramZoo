library(ConsensusPeaks)

# Take a vector of values and get the histogram for integer breaks
obs.to.int.hist = function(x) {
  # Pad with ensured 0 on each side so that we will always have a local minimum
  #   on the end points
  table(cut(x, breaks = (floor(min(x)) - 1):(ceiling(max(x)) + 1)))
}

# Plots the vector x of counts (or table) and the optional segment points s
plot.segments = function(x, s = NULL) {
  index = seq_along(x)
  if(!is.null(s)) {
    opar = par(mfrow = c(2,1), mar = c(2,2,2,2))
  }
  plot(x)

  minmax = local.minmax(x)
  min.ind = minmax$min.ind
  max.ind = minmax$max.ind

  points(seq_along(x)[min.ind], x[min.ind], col = "green")
  points(seq_along(x)[max.ind], x[max.ind], col = "red")

  if(!is.null(s)) {
    plot(x)
    points(s, x[s], col = "orange")
    par(opar)
  }
}

set.seed(13)
# Example 1, uniform
x = obs.to.int.hist(runif(1000, min = 0, max = 50))
s = ftc(x)
plot.segments(x, s)

#set.seed(13)
# set.seed(26)
# Example 2, normal dist
x.norm = obs.to.int.hist(rnorm(1000) * 5)
s.norm = ftc(x.norm)
plot.segments(x.norm, s.norm)

set.seed(26)
# Example 2.1, normal dist
x.norm = obs.to.int.hist(rnorm(1000) * 5)
s.norm = ftc(x.norm)
plot.segments(x.norm, s.norm)

# Example 1, uniform
s = ftc.helen(x)
plot.segments(x, s)

#set.seed(13)
# set.seed(26)
# Example 2, normal dist
x.norm = obs.to.int.hist(rnorm(1000) * 5)
s.norm = ftc.helen(x.norm)
plot.segments(x.norm, s.norm)

set.seed(26)
# Example 2.1, normal dist
x.norm = obs.to.int.hist(rnorm(1000) * 5)
s.norm = ftc.helen(x.norm)
plot.segments(x.norm, s.norm)


# Example 2, messy mixture of gaussians and uniform

set.seed(2623333)
# x.norm.mix = obs.to.int.hist(c(rnorm(25, mean = 1), rnorm(10, mean = 5), c(runif(20, min = 3, max = 5))))
x.norm.mix = obs.to.int.hist(c(rnorm(25, mean = 1), rnorm(80, mean = 8)))
(s.norm.mix.helen = ftc.helen(x.norm.mix))
(s.norm.mix = ftc(x.norm.mix))
plot.segments(x.norm.mix, s.norm.mix.helen)
plot.segments(x.norm.mix, s.norm.mix)
x.norm.mix
