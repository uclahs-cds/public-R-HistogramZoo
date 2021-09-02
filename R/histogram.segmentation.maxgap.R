
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
  c(mint = relative.entropy >= (1/N)*log(L*(L+1)/2) && prob.diff, entropy = relative.entropy)
}

# Meaningful gap
meaningful.gap = function(h, p, a, b, N, L){
  relative.entropy = rel.entropy(h, p, a, b)
  prob.diff = calc.prob.diff(h, p, a, b)
  mgap = (relative.entropy >= (1/N)*log(L*(L+1)/2) && !prob.diff) || (all(h == 0))
  c(mgap = mgap, entropy = relative.entropy)
}

maximal.meaningful = function(x) {
  curr.df = x
  max.intervals = data.frame()
  while(nrow(curr.df) > 0) {
    max.entropy.index = which.max(curr.df$entropy)
    max.entropy = curr.df[max.entropy.index, ]
    max.entropy.seq = seq(max.entropy$Var1, max.entropy$Var2)
    max.intervals = rbind(max.intervals, max.entropy)

    # Find all of the segments that overlap.
    # These will all be less than the maximum
    overlap.max = mapply(function(from, to) {
      s = seq(from, to)
      length(intersect(max.entropy.seq, s)) > 0
    }, from = curr.df$Var1, to = curr.df$Var2)

    curr.df = curr.df[!overlap.max, ]
  }
  # Preserve the old index
  max.intervals$index = rownames(max.intervals)
  rownames(max.intervals) <- NULL
  max.intervals
}

find.all.meaningful.gap = function(x, change.points) {
  # span = seq_along(x)
  # todo = expand.grid(span, span)
  todo = expand.grid(change.points, change.points)
  # todo$Var2 = todo$Var2 - 1 # If there are segments
  todo = todo[todo$Var2 > todo$Var1,]

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
  df = cbind(todo, mgap)
  # df = df[order(df$Var2, df$Var1),]
  df = df[order(df$entropy), ]
  df$mgap = as.numeric(df$mgap)
  # df$scaled_entropy = (df$entropy - min(df$entropy, na.rm = T)) / (max(df$entropy, na.rm = T) - min(df$entropy, na.rm = T))

  seg.gap.data = df[df$mgap > 0 & !is.na(df$mgap), ]
  # seg.gap.data$index = 1:nrow(seg.gap.data)

  maximal.meaningful(seg.gap.data)
}

# Finds the meaningful gaps between the points in s
meaningful.gaps.local = function(x, seg.points, change.points) {

  max.gaps.list <- lapply(seq(2, length(seg.points)), function(i) {
    x.sub = x[seg.points[i-1]:seg.points[i]]
    chg.pts = change.points[change.points >= seg.points[i-1] & change.points <= seg.points[i]] - seg.points[i-1] + 1

    max.gaps = find.all.meaningful.gap(x.sub, chg.pts)

    if(nrow(max.gaps) > 0) {
      max.gaps[, c('Var1','Var2')] <- max.gaps[, c('Var1','Var2')] + seg.points[i-1] -1
      max.gaps$seg.start = seg.points[i-1]
      max.gaps$seg.end = seg.points[i]
      max.gaps
    }
  })

  max.gaps <- do.call(rbind.data.frame, max.gaps.list)
  max.gaps
}

find.new.segments = function(gaps.df) {
  new.segments = c()
  for(i in rownames(gaps.df)) {
    r = gaps.df[i, ]
    # Have a large gap
    if(r$Var2 - r$Var1 > 1) {
      left.diff = abs(r$seg.start - r$Var1)
      right.diff = abs(r$seg.end-r$Var2)
      if(left.diff < right.diff && right.diff > 2) {
        new.segments = c(new.segments, r$Var2 - 1)
      } else if(left.diff > right.diff && left.diff > 2){
        new.segments = c(new.segments, r$Var1 + 1)
      }
    }
  }
  new.segments
}

maximal.gaps.multiplot = function(x, seg.points, max.gaps, new.segments = NULL, ...) {
  # Add the gap data to a vector for heatmap
  maximal.data <- rep("none", length(x))
  for(i in rownames(max.gaps)) {
    r = max.gaps[i,]
    gap.seq = seq(r$Var1, r$Var2) - 1
    maximal.data[gap.seq] <- "gap"
  }
  maximal.data <- factor(maximal.data, levels = c("none", "gap"))

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

  colour.scheme = c('white', 'orange')
  maximal.heatmap <- create.heatmap(
    x = t(as.matrix(as.numeric(maximal.data))),
    clustering.method = 'none',
    scale.data = FALSE,
    colour.scheme = colour.scheme,
    total.colours = 1 + length(colour.scheme),
    grid.col = TRUE,
    print.colour.key = FALSE,
    force.grid.col = TRUE,
    axes.lwd = 0
  )

  segment.data <- rep("none", length(x))
  segment.data[seg.points] <- "segment"
  if(is.null(new.segments)) {
    colour.scheme.segments <- c('white', 'red')
    segment.data <- factor(segment.data, levels = c("none", "segment"))
  } else {
    colour.scheme.segments <- c('white', 'red', 'blue')
    segment.data[new.segments] <- "extra_segment"
    segment.data <- factor(segment.data, levels = c("none", "segment", "extra_segment"))
  }

  segments.heatmap <- create.heatmap(
    x = t(as.matrix(as.numeric(segment.data))),
    clustering.method = 'none',
    scale.data = FALSE,
    colour.scheme = colour.scheme.segments,
    total.colours = 1 + length(colour.scheme.segments),
    grid.col = TRUE,
    print.colour.key = FALSE,
    force.grid.col = TRUE,
    axes.lwd = 0
  )

  create.multipanelplot(
    plot.objects = list(bp, segments.heatmap, maximal.heatmap),
    plot.objects.heights = c(1, 0.15, 0.15),
    y.spacing = -0.75,
    ylab.label = "Segments",
    xlab.label = "Index",
    ylab.cex = 2,
    xlab.cex = 2,
    height = 6,
    width = 8,
    ...
  )
}
