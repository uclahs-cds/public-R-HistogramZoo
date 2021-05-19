library(ConsensusPeaks)

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
set.seed(26)
# x.norm.mix = obs.to.int.hist(c(rnorm(25, mean = 1), rnorm(10, mean = 5), c(runif(20, min = 3, max = 5))))
x.norm.mix = obs.to.int.hist(c(rnorm(25, mean = 1), rnorm(80, mean = 10)))
s.norm.mix.helen = ftc.helen(x.norm.mix)
s.norm.mix = ftc(x.norm.mix)
plot.segments(x.norm.mix, s.norm.mix.helen)
plot.segments(x.norm.mix, s.norm.mix)
x.norm.mix
