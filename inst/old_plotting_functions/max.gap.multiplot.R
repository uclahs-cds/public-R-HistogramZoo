
maximal.gaps.multiplot = function(x, seg.points, max.gaps, new.segments = NULL, ...) {
  # Add the gap data to a vector for heatmap
  maximal.data <- rep("none", length(x))
  for(i in rownames(max.gaps)) {
    r = max.gaps[i,]
    gap.seq = seq(r$start, r$end) - 1
    maximal.data[gap.seq] <- "gap"
  }
  maximal.data <- factor(maximal.data, levels = c("none", "gap"))

  bp = create.barplot(
    Freq ~ start,
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
