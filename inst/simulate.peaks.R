set.seed(26)
x.norm.mix = obs.to.int.hist(c(rnorm(25, mean = 1), rnorm(80, mean = 10), runif(50, min = 25, max = 35)))

res = segment.fit.agnostic(x.norm.mix)

res.trunc = segment.fit.agnostic(x.norm.mix, truncated.models = TRUE, remove.low.entropy =T)

new.s = unique(unlist(res[res$final == 1, c('seg.start', 'seg.end')], use.names = F))
new.s.bin = names(x.norm.mix)[new.s]
plot.segments(x.norm.mix, s = new.s.bin)

res[res$final == 1, ]
res.trunc[res.trunc$final == 1, ]
