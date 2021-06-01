library(ConsensusPeaks)

# I don't think this is a bug in the code, but rather a strange result
# Two clear distinct segments that are uniform but gets merged into one segment
set.seed(123)
x.long.unif <- obs.to.int.hist(runif(100, max = 50))
x.long.unif[1:(floor(length(x.long.unif) /  2))] <- x.long.unif[1:(floor(length(x.long.unif) /  2))] + 10
s.long.unif <- ftc(x.long.unif)
s.long.unif.helen <- ftc.helen(x.long.unif, s = local.minmax(x.long.unif))
plot.segments(x.long.unif, s.long.unif)
# plot.segments(x.long.unif, s.long.unif.helen)


# Failing example for non min/max
x <- c(10,10,0,0,0,0,10,10)
names(x) <- seq_along(x)
s <- ftc.helen(x, s = local.minmax(x))
plot.segments(x, unlist(s))
s2 <- ftc(x)

plot.segments(x, s2)


# Example showing why we "need" to merge J > 1
set.seed(21)
x <- obs.to.int.hist(unlist(lapply(seq(1, 100, by = 4), function(x) rnorm(100, mean = x))))
s <- ftc.helen(x, s = local.minmax(x))
s2 <- ftc(x)
plot.segments(x, unlist(s))
plot.segments(x, s2)

# Failing example for min/max method
x.fail <- structure(c(`(-1,0]` = 0L, `(0,1]` = 16L, `(1,2]` = 16L, `(2,3]` = 6L,
                      `(3,4]` = 3L, `(4,5]` = 0L, `(5,6]` = 4L, `(6,7]` = 1L, `(7,8]` = 1L,
                      `(8,9]` = 1L, `(9,10]` = 1L, `(10,11]` = 1L, `(11,12]` = 0L), .Dim = 13L, .Dimnames = structure(list(
                        c("(-1,0]", "(0,1]", "(1,2]", "(2,3]", "(3,4]", "(4,5]",
                          "(5,6]", "(6,7]", "(7,8]", "(8,9]", "(9,10]", "(10,11]",
                          "(11,12]")), .Names = ""), class = "table")
s.fail.helen = ftc.helen(x.fail, s = local.minmax(x.fail))
# This fails for min/max method because the min/maxes are not alternated correctly
s.fail = ftc(x.fail)
plot.segments(x.fail, unlist(s.fail.helen))
plot.segments(x.fail, unlist(s.fail))
