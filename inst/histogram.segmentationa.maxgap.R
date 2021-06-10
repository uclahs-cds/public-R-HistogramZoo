
library(isotone)
library(ggplot2)
library(BoutrosLab.plotting.general)
library(reshape2)

# Preamble ----------------------------------------------------------------
# Trying to implement a gap detection thing based on the Delon papers


# Functions ---------------------------------------------------------------

# Take a vector of values and get the histogram for integer breaks
obs.to.int.hist = function(x) {
  table(cut(x, breaks = (floor(min(x)) - 1):(ceiling(max(x)) + 1)))
}

# Relative entropy
rel.entropy = function(h, p, a, b) {
  interval = a:b
  # Round to prevent floating point issues
  hab = round(sum(h[interval]), digits = 14)
  pab = round(sum(p[interval]), digits = 14)
  if(pab == 0 || pab == 1 || hab == 0 || hab == 1) {
    return(0)
  } else {
    hab * log(hab / pab) + (1 - hab) * log((1 - hab) / (1 - pab))
  }
}

# Plots the vector x of counts (or table) and the optional segment points s
plot.segments = function(x, s = NULL, threshold = 0, ...) {
  index = seq_along(x)
  plot(x, type = "h", ...)
  
  if(!is.null(s)) {
    points(s, x[s], col = "orange")
  }
}

# Uniform Generation
generate.unif = function(x){
  rep(1/length(x), length(x))
}

calc.prob.diff = function(h, p, a, b){
  interval = a:b
  # Round to prevent floating point issues
  hab = round(sum(h[interval]), digits = 14)
  pab = round(sum(p[interval]), digits = 14)
  hab > pab
}

# Meaningful interval
meaningful.interval = function(h, p, a, b, N, L){
  relative.entropy = rel.entropy(h, p, a, b)
  prob.diff = calc.prob.diff(h, p, a, b)
  relative.entropy >= (1/N)*log(L*(L+1)/2) && prob.diff
}

# Meaningful gap
meaningful.gap = function(h, p, a, b, N, L){
  relative.entropy = rel.entropy(h, p, a, b)
  prob.diff = calc.prob.diff(h, p, a, b)
  relative.entropy >= (1/N)*log(L*(L+1)/2) && !prob.diff
}

# Test Cases --------------------------------------------------------------

set.seed(5)
n1 = rnorm(100, mean = 5, sd = 2)
n2 = rnorm(50, mean = 10, sd = 3)
n3 = rnorm(200, mean = 30, sd = 1)
n.combined = c(n1, n2, n3)
x = obs.to.int.hist(n.combined)
plot.segments(x)

# Running FTC
s <- ftc.helen(x, s = sort(unlist(local.minmax(x))))
s2 <- ftc(x)
plot.segments(x, unlist(s))
plot.segments(x, s2)

# Trying Maximal Meaningful
todo = expand.grid(c(1:35), c(1:35))
todo = todo[todo$Var2 > todo$Var1,]
mint = do.call(rbind, lapply(1:nrow(todo), function(i) {
  meaningful.interval(
    h = x/sum(x),
    p = generate.unif(x),
    a = todo$Var1[i],
    b = todo$Var2[i],
    N = sum(x),
    L = length(x)
  )
}))  
mgap = do.call(rbind, lapply(1:nrow(todo), function(i) {
  meaningful.gap(
    h = x/sum(x),
    p = generate.unif(x),
    a = todo$Var1[i],
    b = todo$Var2[i],
    N = sum(x),
    L = length(x)
  )
})) 
df = cbind(todo, mint, mgap)
df = df[order(df$Var1, df$Var2),]
df$mint = as.numeric(df$mint)
df$mgap = as.numeric(df$mgap)

# Create A Heatmap Looking For Meaningful Segments
plotting.mint = acast(df, Var1 ~ Var2, value.var = "mint")
plotting.mgap = acast(df, Var1 ~ Var2, value.var = "mgap")

bp = create.barplot(
  Freq ~ Var1,
  data.frame(x),
  xaxis.cex = 0,
  xlab.cex = 0,
  ylab.label = "Counts",
  xaxis.tck = 0,
  yaxis.tck = 0,
  yaxis.cex = 0.8,
  ylab.cex = 1
)

hm.int = create.heatmap(
  plotting.mint,
  clustering.method = "none",
  ylab.label = "Meaningful Intervals",
  xaxis.cex = 0.8,
  yaxis.cex = 0.8,
  ylab.cex = 1,
  xaxis.lab = colnames(plotting.mint),
  yaxis.lab = rownames(plotting.mint),
  xaxis.rot = F,
  colour.scheme = c("white", "black"),
  same.as.matrix = T,
  print.colour.key = F
)

hm.gap = create.heatmap(
  plotting.mgap,
  clustering.method = "none",
  ylab.label = "Meaningful Gaps",
  xaxis.cex = 0.8,
  yaxis.cex = 0.8,
  ylab.cex = 1,
  xaxis.lab = colnames(plotting.mgap),
  yaxis.lab = rownames(plotting.mgap),
  xaxis.rot = F,
  colour.scheme = c("white", "black"),
  same.as.matrix = T,
  print.colour.key = F
)

filename = "~/Desktop/Meaningful.Segments.pdf"
pdf(filename, height = 12, width = 8)
create.multipanelplot(
  plot.objects = list(bp, hm.int, hm.gap),
  plot.objects.heights = c(1, 4, 4),
  y.spacing = -1,
  ylab.label = "Start",
  xlab.label = "End",
  ylab.cex = 2,
  xlab.cex = 2
)
dev.off()