
bpg.plot = function(
  output.tag,
  output.dir,
  distr.plotting.data,
  geneinfo,
  bin.counts,
  p,
  results
) {

  # metric.label = switch(
  #   met,
  #   "jaccard" = "Jaccard Index",
  #   "intersection" = "Histogram Intersection",
  #   "ks" = "KS Statistic",
  #   "mse" = "Mean Squared Error",
  #   "chisq" = "Chi-square statistic")

  # Bin counts
  bin.counts = bin.counts[,c("start", "Coverage")]
  bin.counts$dist = "coverage"
  colnames(bin.counts) = c('x', 'dens', 'dist')

  # Plotting data
  distr.plotting.data = c(distr.plotting.data, list(bin.counts))
  for(i in 1:length(distr.plotting.data)){distr.plotting.data[[i]]$'i' <- i}
  lineplot.data = do.call(rbind, distr.plotting.data)

  # Points
  points = bin.counts[bin.counts$x %in% p,]

  # Residuals
  fitted.distribution.density = lineplot.data[lineplot.data$dist != "coverage",]
  residual.data = merge(fitted.distribution.density[,c("x", "dens")], bin.counts[,c("x", "dens")], by = "x", all = T)
  colnames(residual.data) = c("x", "fitted", "real")
  residual.data$resid = residual.data$real - residual.data$fitted

  # Reference Colours, Numeric & Line Width for Distributions
  col.reference = structure(c("black", "darkorange", "chartreuse4", "darkorchid4"), names = c("coverage", "norm", "gamma", "unif"))
  col.numeric = structure(1:3, names = c("norm", "gamma", "unif"))
  lwd.reference = structure(c(1, 2.5, 2.5, 2.5), names = c("coverage", "norm", "gamma", "unif"))

  # Plotting attributes for lineplot
  col.used = unique(lineplot.data[,c("dist", "i")])
  col.used = col.used[order(col.used$i),]
  col.vec = col.reference[col.used$dist]
  lwd.vec = lwd.reference[col.used$dist]

  # Plotting
  sc = BoutrosLab.plotting.general::create.scatterplot(
    dens ~ x,
    lineplot.data,
    # Groups
    groups = lineplot.data$i,
    col = col.vec,
    # Axes
    xlimits = c(0, geneinfo$exome_length),
    xaxis.cex = 0,
    xlab.cex = 0,
    ylab.cex = 1,
    xaxis.tck = 0,
    yaxis.cex = 0.8,
    ylab.label = "Coverage (at BP resolution)",
    main.cex = 0,
    # Lines & PCH
    type = c('a'),
    lwd = lwd.vec,
    # Adding extra points
    add.points = T,
    points.x = points$x,
    points.y = points$dens,
    points.pch = 19,
    points.col = 'red'
  )

  resid.chgpts = which(abs(diff(sign(residual.data$resid))) == 2)
  sc.res =  BoutrosLab.plotting.general::create.scatterplot(
    resid ~ x,
    residual.data,
    # Groups
    # groups = residual.data$something,
    col = "black",
    # Axes
    xlimits = c(0, geneinfo$exome_length),
    xaxis.cex = 0,
    xlab.cex = 0,
    ylab.cex = 1,
    ylab.label = "Residuals",
    xaxis.tck = 0,
    yaxis.cex = 0.8,
    main.cex = 0,
    # Lines & PCH
    type = c('p'),
    cex = 0.01,
    # Adding lines at changepoints
    abline.v = resid.chgpts,
    abline.h = 0,
    abline.lty = "dotted",
    abline.col = "lightgrey",
    abline.lwd = 0.01
  )

  # Goodness of Fit
  gof.mat = matrix(NA, ncol = 1, nrow = geneinfo$exome_length, dimnames = list(1:geneinfo$exome_length, "val"))
  res.final = results[results$final == 1,, drop = F]
  for(i in 1:nrow(res.final)){gof.mat[res.final$seg.start[i]:res.final$seg.end[i],1] <- res.final$value[i]}
  # Legend metrics
  min.at = min(min(gof.mat, na.rm = T), 0.8)
  max.at = 1

  gof.hm = BoutrosLab.plotting.general::create.heatmap(
      gof.mat[,c("val", "val")],
      clustering.method = 'none',
      # Plotting Characteristics
      axes.lwd = 0,
      yaxis.tck = 0,
      # Y axis Labels
      yaxis.cex = 0.8,
      yaxis.lab = "Jaccard Index",
      yat = 1.5,
      # Colours
      colour.scheme = c("dodgerblue4", "gold"),
      at = seq(min.at, max.at, length.out = 10),
      fill.colour = "white",
      # Adding lines for segments
      force.grid.col = TRUE,
      grid.col = TRUE,
      col.lines = peak.endpoints,
      col.lwd = 1,
      print.colour.key = F
  )

  # Distribution Voting
  peak.endpoints = sort(unique(c(results$seg.start, results$seg.end)))
  results$metric[results$final == 1] <- "MaxVote"
  unique.mets = unique(results$metric)
  unique.mets = unique.mets[order(match(unique.mets, rev(c("MaxVote", "jaccard", "intersection", "mse", "chisq", "ks"))))]
  result.mat = matrix(NA, ncol = length(unique.mets), nrow = geneinfo$exome_length, dimnames = list(1:geneinfo$exome_length, unique.mets))
  for(i in 1:nrow(results)){result.mat[results$seg.start[i]:results$seg.end[i],results$metric[i]] <- col.numeric[results$dist[i]]}

  dist.hm = BoutrosLab.plotting.general::create.heatmap(
    result.mat,
    clustering.method = 'none',
    # Plotting Characteristics
    axes.lwd = 0,
    yaxis.tck = 0,
    # Y axis Labels
    yaxis.cex = 0.8,
    yaxis.lab = colnames(result.mat),
    ylab.cex = 1,
    ylab.label = "Metrics",
    # Colours
    at = seq(-0.5, length(col.reference), 1),
    total.colours = length(col.reference) + 1,
    colour.scheme = c('white', col.reference[2:length(col.reference)]),
    fill.colour = "white",
    # Adding lines for segments
    force.grid.col = TRUE,
    grid.col = TRUE,
    col.lines = peak.endpoints,
    col.lwd = 1,
    force.grid.row = TRUE,
    grid.row = T,
    row.lines = 1:ncol(result.mat)+0.5,
    row.lwd = 1,
    row.colour = c(rep("black", ncol(result.mat)-2), "red"),
    # Colourkey
    print.colour.key = F
  )

  # Generating Exons
  genegr = GenomicRanges::GRanges(seqnames = geneinfo$chr, IRanges::IRanges(start = 1, end = geneinfo$exome_length), strand = geneinfo$strand)
  anno = geneinfo$anno
  anno$start = geneinfo$DNA2RNA[anno$start - geneinfo$left+1]
  anno$stop = geneinfo$DNA2RNA[anno$stop - geneinfo$left+1]
  anno.gr = GenomicRanges::makeGRangesFromDataFrame(anno, keep.extra.columns = T)
  anno.gr = S4Vectors::split(anno.gr, anno.gr$transcript)
  anno.coverage = lapply(anno.gr, GenomicRanges::coverage)
  bins = unlist(GenomicRanges::tile(genegr, width = 1))
  anno.counts = lapply(anno.coverage, function(i) GenomicRanges::binnedAverage(bins, i, "Coverage"))
  transcript.coverage = do.call(cbind.data.frame, lapply(anno.counts, function(x) x$Coverage))

  x.digits = floor(log10(geneinfo$exome_length))-1
  x.at = seq(0, round(geneinfo$exome_length, digits = -x.digits), length.out = 5)
  hm.coverage = BoutrosLab.plotting.general::create.heatmap(
    transcript.coverage,
    clustering.method = 'none',
    # Axes labels
    xaxis.lab = x.at,
    xat = x.at,
    xaxis.rot = 0,
    xaxis.cex = 1,
    xaxis.tck = 1,
    yaxis.tck = 0,
    xlab.cex = 1,
    xlab.label = "Transcript Coordinate",
    # Discrete Colours
    at = c(-0.5, 0.5, 1.5),
    total.colours = 3,
    colour.scheme = c('white', 'pink'),
    # Border
    axes.lwd = 1,
    # Adding lines for segments
    force.grid.col = TRUE,
    grid.col = TRUE,
    col.lines = peak.endpoints,
    col.lwd = 1,
    force.grid.row = TRUE,
    grid.row = TRUE,
    row.lines = 1:ncol(transcript.coverage)+0.5,
    row.lwd = 1,
    # Colourkey
    print.colour.key = F
  )

  # Legend
  covariate.legend <- list(
    legend = list(
      colours = 'red',
      labels = c("FTC Endpoints"),
      title = expression(bold(underline('Points'))),
      lwd = 0.5
    ),
    legend = list(
      colours = 'black',
      labels = 'Peak Coverage',
      title = expression(bold(underline('Lines'))),
      lwd = 0.5
    ),
    legend = list(
      colours = c('black', 'lightgrey'),
      labels = c('Obs - Exp', 'Intersect points'),
      title = expression(bold(underline('Residuals'))),
      lwd = 0.5
    ),
    legend = list(
      colours = col.reference[2:length(col.reference)],
      labels = names(col.reference)[2:length(col.reference)],
      title = expression(bold(underline('Distributions'))),
      lwd = 0.5
    ),
    legend = list(
      colours = c('dodgerblue4', 'gold'),
      labels = c(formatC(min.at, format = "g", digits = 2), max.at),
      title = bquote(bold(underline(.("Jaccard Index")))),
      continuous = TRUE,
      height = 2,
      angle = -90,
      tck = 1,
      tck.number = 3,
      at = c(0,100),
      labels.rot = 0,
      lwd = 0.5
    ),
    legend = list(
      colours = c("pink"),
      labels = c("Exon"),
      title = expression(bold(underline('Transcript Coverage'))),
      lwd = 0.5
    )
  )

  side.legend <- BoutrosLab.plotting.general::legend.grob(
    legends = covariate.legend,
    label.cex = 0.7,
    title.cex = 0.7,
    title.just = 'left',
    title.fontface = 'bold',
    size = 2
  )

  # Plotting
  filename = file.path(output.dir, paste0(geneinfo$gene, ".", output.tag, ".SegmentAndFit.pdf"))
  pdf(filename, width = 10, height = 10)

  transcript.height = min(3, ncol(transcript.coverage)*0.5) + 0.55
  distplot.height = min(2.8, length(unique.mets)*0.8)

  mpp = BoutrosLab.plotting.general::create.multipanelplot(
    plot.objects = list(sc, sc.res, gof.hm, dist.hm, hm.coverage),
    plot.objects.heights = c(10, 3, 1, distplot.height, transcript.height),
    y.spacing = -5,
    # Labels
    main = geneinfo$gene,
    main.cex = 2,
    # Legend
    legend = list(
      right = list(
        x = 0.8,
        y = 1,
        fun = side.legend
      )
    )
  )

  print(mpp)
  dev.off()

}
