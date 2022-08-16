
library(HistogramZoo)


# Preamble ----------------------------------------------------------------
# Thinking of switching find_stepfunction_changepoints to find_local_optima
# in segment_and fit


# Testing for Stepfunctions -----------------------------------------------

# Generating data
data = c(rep(1, 3), rep(2, 4), rep(3, 3), rep(2, 1), rep(1, 5))
xhist = Histogram(histogram_data = data)
create_coverageplot(xhist)

# changepts - this finds the lower of the points
chgpts = HistogramZoo:::find_stepfunction_chgpts(data)
create_coverageplot(xhist, add.points = T, points.x = chgpts, points.y = data[chgpts], points.col = "red")

# find_local_optima
optima = find_local_optima(data)
optima = sort(unlist(optima))
create_coverageplot(xhist, add.points = T, points.x = optima, points.y = data[optima], points.col = "red")

# Testing for Smooth data -------------------------------------------------

set.seed(314)
dt = rnorm(10000, mean = 20, sd = 10)
xhist = observations_to_histogram(dt, histogram_bin_width = 2)
interval_start = xhist$interval_start
data = xhist$histogram_data
create_coverageplot(xhist)

# changepts - this finds the lower of the points
chgpts = HistogramZoo:::find_stepfunction_chgpts(data)
create_coverageplot(xhist, add.points = T, points.x = interval_start[chgpts], points.y = data[chgpts], points.col = "red")

# find_local_optima
optima = find_local_optima(data)
optima = sort(unlist(optima))
create_coverageplot(xhist, add.points = T, points.x = interval_start[optima], points.y = data[optima], points.col = "red")

