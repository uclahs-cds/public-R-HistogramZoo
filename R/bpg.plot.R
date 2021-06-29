
bpg.plot = function(
  output.tag,
  output.dir,
  distr.plotting.data,
  geneinfo,
  bin.counts,
  seg.gr,
  p,
  fitted.seg.gr,
  results
) {

  # Plotting data
  n.distributions = length(distr.plotting.data)
  for(i in 1:n.distributions){distr.plotting.data[[i]]$'i' <- i}
  distplot.data = do.call(rbind, distr.plotting.data)

  # Bin counts
  bin.counts = bin.counts[,c("start", "Coverage")]
  bin.counts$dist = "coverage"
  bin.counts$i = n.distributions + 1
  colnames(bin.counts) = c('x', 'dens', 'dist', 'i')

  # Points
  points = bin.counts[bin.counts$x %in% p,]

  # Plotting Data
  lineplot.data = rbind(distplot.data, bin.counts)

  # Colours
  col.reference = structure(c("black", "darkorange", "chartreuse4", "darkorchid4"),
                            names = c("coverage", "norm", "gamma", "unif"))
  col.numeric = structure(1:3, names = c("norm", "gamma", "unif"))
  lwd.reference = structure(c(1, 2.5, 2.5, 2.5), names = c("coverage", "norm", "gamma", "unif"))
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
    yaxis.cex = 1,
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

  # Generating Gene Segments
  genegr = GenomicRanges::GRanges(seqnames = geneinfo$chr, IRanges::IRanges(start = 1, end = geneinfo$exome_length), strand = geneinfo$strand)
  anno = geneinfo$anno
  anno$start = geneinfo$DNA2RNA[anno$start - geneinfo$left+1]
  anno$stop = geneinfo$DNA2RNA[anno$stop - geneinfo$left+1]
  anno.gr = GenomicRanges::makeGRangesFromDataFrame(anno, keep.extra.columns = T)
  anno.gr = S4Vectors::split(anno.gr, anno.gr$transcript)
  anno.coverage = lapply(anno.gr, GenomicRanges::coverage)
  bins = unlist(GenomicRanges::tile(genegr, width = 1))
  anno.counts = lapply(anno.coverage, function(i) GenomicRanges::binnedAverage(bins, i, "Coverage"))
  transcript.coverage = do.call(cbind, lapply(anno.counts, function(x) x$Coverage))
  transcript.coverage = as.data.frame(transcript.coverage, stringsAsFactors = F)

  # Adding p-values & segments
  # Change seg.gr to fitted.seg.gr for just the fitted segments
  ovl = GenomicRanges::findOverlaps(bins, seg.gr)
  bins$metric = NA
  bins$metric[S4Vectors::queryHits(ovl)] = results$metric[S4Vectors::subjectHits(ovl)]

  # Adding segment fits
  ovl = GenomicRanges::findOverlaps(bins, fitted.seg.gr)
  bins$fitted = 0
  bins$fitted[S4Vectors::queryHits(ovl)] = col.numeric[results$dist[S4Vectors::subjectHits(ovl)]]

  # Plotting Heatmap
  heatmap.data = data.frame(bins, stringsAsFactors = F)
  heatmap.data = heatmap.data[,c("start", "fitted", "metric")]
  heatmap.data = heatmap.data[order(heatmap.data$start),]

  hm.fu = BoutrosLab.plotting.general::create.heatmap(
    heatmap.data[,c("fitted", "fitted")],
    clustering.method = 'none',
    # Plotting Characteristics
    axes.lwd = 1,
    yaxis.tck = 0,
    # Discrete Colours
    at = seq(-0.5, length(col.reference), 1),
    total.colours = length(col.reference) + 1,
    colour.scheme = c('white', col.reference[2:length(col.reference)]),
    # Adding lines for segments
    force.grid.col = TRUE,
    grid.col = TRUE,
    col.lines = c(GenomicRanges::start(fitted.seg.gr), GenomicRanges::end(fitted.seg.gr)),
    # Colourkey
    # colourkey.labels.at = seq(0, length(seg.gr)+1, 1),
    # colourkey.labels = seq(0, length(seg.gr)+1, 1),
    print.colour.key = F
  )

  hm.gof = BoutrosLab.plotting.general::create.heatmap(
    heatmap.data[,c("metric", "metric")],
    clustering.method = 'none',
    # Plotting Characteristics
    axes.lwd = 1,
    yaxis.tck = 0,
    # Colours
    colour.scheme = c('dodgerblue4', 'firebrick'),
    at = seq(0, 1, 0.05), # fix this when scaled
    fill.colour = "white",
    # Adding lines for segments
    force.grid.col = TRUE,
    grid.col = TRUE,
    col.lines = c(GenomicRanges::start(seg.gr), GenomicRanges::end(seg.gr)),
    col.lwd = 1,
    # Colourkey
    print.colour.key = F
  )

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
    col.lines = c(GenomicRanges::start(seg.gr), GenomicRanges::end(seg.gr)),
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
      colours = col.reference[2:length(col.reference)],
      labels = names(col.reference)[2:length(col.reference)],
      title = expression(bold(underline('Distributions'))),
      lwd = 0.5
    ),
    legend = list(
      colours = c('dodgerblue4', 'firebrick'),
      labels = c(0, 1),
      title = expression(bold(underline('Jaccard Index'))),
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

  mpp = BoutrosLab.plotting.general::create.multipanelplot(
    plot.objects = list(sc, hm.fu, hm.gof, hm.coverage),
    plot.objects.heights = c(10, 1, 1, transcript.height),
    y.spacing = -4,
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
