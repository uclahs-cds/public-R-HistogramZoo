
coverage.plot = function(
  histogram.coverage,
  dist.data,
  p
){
  # Plotting data
  xlim = nrow(histogram.coverage)
  points = histogram.coverage[histogram.coverage$x %in% p,]
  plt.data = c(dist.data, list(histogram.coverage))
  for(i in 1:length(plt.data)){plt.data[[i]]$'i' <- i}
  plt.data = do.call(rbind, plt.data)

  # Reference Colours, Numeric & Line Width for Distributions
  col.reference = c("black", "darkorange", "chartreuse4", "darkorchid4")
  names(col.reference) = c("coverage", "norm", "gamma", "unif")
  lwd.reference = c(1, 2.5, 2.5, 2.5)
  names(lwd.reference) = c("coverage", "norm", "gamma", "unif")

  # Plotting attributes for lineplot
  col.used = unique(plt.data[,c("dist", "i")])
  col.used = col.used[order(col.used$i),]
  col.vec = col.reference[col.used$dist]
  lwd.vec = lwd.reference[col.used$dist]

  plt = BoutrosLab.plotting.general::create.scatterplot(
    dens ~ x,
    plt.data,
    # Groups
    groups = plt.data$i,
    col = col.vec,
    # Axes
    xlimits = c(0, xlim),
    xaxis.cex = 0,
    xlab.cex = 0,
    ylab.cex = 1,
    xaxis.tck = 0,
    yaxis.cex = 0.8,
    ylab.label = "Coverage",
    main.cex = 0,
    yaxis.fontface = 1,
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

  return(plt)
}

residual.plot = function(
  histogram.coverage,
  dist.data
){

  xlim = nrow(histogram.coverage)
  dist.data = do.call(rbind, dist.data)
  residual.data = merge(dist.data[,c("x", "dens")], histogram.coverage[,c("x", "dens")], by = "x", all = T)
  colnames(residual.data) = c("x", "fitted", "real")
  residual.data$resid = residual.data$real - residual.data$fitted

  residual.chgpts = which(abs(diff(sign(residual.data$resid))) == 2)
  plt =  BoutrosLab.plotting.general::create.scatterplot(
    resid ~ x,
    residual.data,
    # Groups
    # groups = residual.data$something,
    col = "black",
    # Axes
    xlimits = c(0, xlim),
    xaxis.cex = 0,
    xlab.cex = 0,
    ylab.cex = 1,
    ylab.label = "Residuals",
    xaxis.tck = 0,
    yaxis.cex = 0.8,
    main.cex = 0,
    yaxis.fontface = 1,
    # Lines & PCH
    type = c('p'),
    cex = 0.01,
    # Adding lines at changepoints
    abline.v = residual.chgpts,
    abline.h = 0,
    abline.lty = "dotted",
    abline.col = "lightgrey",
    abline.lwd = 0.01
  )
  return(plt)
}

metrics.plot = function(
  models,
  x.limit
){

  # Goodness of Fit
  init = rep(NA, x.limit)
  for(i in seq_along(models)){
    m = models[[i]][['majority.vote']]
    init[m$seg.start:m$seg.end] = m$value
  }
  final.mat = matrix(init, ncol = 1)

  # Legend metrics
  min.at = min(final.mat, na.rm = T)
  min.at = min(min.at, 0.8)
  max.at = 1

  plt = BoutrosLab.plotting.general::create.heatmap(
    final.mat[,c(1, 1)],
    clustering.method = 'none',
    # Plotting Characteristics
    axes.lwd = 0,
    yaxis.tck = 0,
    # Y axis Labels
    yaxis.cex = 0.8,
    yaxis.lab = "Jaccard Index",
    yat = 1.5,
    yaxis.fontface = 1,
    # Colours
    colour.scheme = c("dodgerblue4", "gold"),
    at = seq(min.at, max.at, length.out = 10),
    fill.colour = "white",
    # Adding lines for segments
    # force.grid.col = TRUE,
    # grid.col = TRUE,
    # col.lines = peak.endpoints,
    # col.lwd = 1,
    print.colour.key = F
  )

  # Colors to numeric
  col.num = 1:3
  names(col.num) =  c("norm", "gamma", "unif")

  # Distribution Voting
  ref.mets = c("Consensus", "jaccard", "intersection", "mse", "chisq", "ks")
  unique.mets = sapply(models[[1]], `[[`, 'metric')
  unique.mets = unique.mets[order(match(unique.mets, rev(ref.mets)))]
  metric.mat = matrix(
    NA,
    ncol = length(unique.mets),
    nrow = x.limit,
    dimnames = list(NULL, unique.mets)
    )
  models = unlist(models, recursive = F)
  for(i in seq_along(models)){
    m = models[[i]]
    metric.mat[m$seg.start:m$seg.end, m$metric] <- col.num[m$dist]
  }

  # Labels
  metric.label = c(expression(bold("Consensus")), "Jaccard", "Intersection", "K-S", "MSE", expression(paste(chi^"2")))
  names(metric.label) = c("Consensus", "jaccard", "intersection", "ks", "mse", "chisq")

  # Colours
  col.reference = c("black", "darkorange", "chartreuse4", "darkorchid4")
  names(col.reference) = c("coverage", "norm", "gamma", "unif")

  # X axis
  x.digits = floor(log10(x.limit))-1
  x.at = seq(0, round(x.limit, digits = -x.digits), length.out = 5)
  x.at[x.at == 0] <- 1

  plt2 = BoutrosLab.plotting.general::create.heatmap(
    metric.mat,
    clustering.method = 'none',
    # Plotting characteristics
    axes.lwd = 0,
    # Y axis labels
    yaxis.lab = metric.label[colnames(metric.mat)],
    yaxis.cex = 0.8,
    ylab.cex = 1,
    yaxis.tck = 0,
    ylab.label = "Metrics",
    # X axis labels
    xaxis.lab = x.at,
    xat = x.at,
    xaxis.rot = 0,
    xaxis.cex = 1,
    xaxis.tck = 1,
    xlab.cex = 1,
    xlab.label = "Transcript Coordinate",
    xaxis.fontface = 1,
    # Colours
    at = seq(-0.5, length(col.reference), 1),
    total.colours = length(col.reference) + 1,
    colour.scheme = c('white', col.reference[2:length(col.reference)]),
    fill.colour = "white",
    # Adding lines for segments
    # force.grid.col = TRUE,
    # grid.col = TRUE,
    # col.lines = peak.endpoints,
    # col.lwd = 1,
    force.grid.row = TRUE,
    grid.row = T,
    row.lines = 1:ncol(metric.mat)+0.5,
    row.lwd = 1,
    row.colour = c(rep("black", ncol(metric.mat)-2), "red"),
    # Colourkey
    print.colour.key = F
  )
  return(list(plt, plt2, min.at))
}

plot.fitted.segments = function(
  histogram.names,
  coverage.model.obj,
  file.name,
  plot.types = c("coverage", "residuals", "metrics", "transcript")
){

  models = lapply(coverage.model.obj$results, `[[`, "models")
  point.vals = lapply(coverage.model.obj$results, `[[`, "p")

  compiled.plts = list()
  # Loop through genes to plot things
  for(i in histogram.names){

    # Histogram coverage
    histogram.coverage = coverage.model.obj$histogram.coverage[[i]]
    xlim = length(histogram.coverage)
    histogram.coverage = data.frame(
      "x" = 1:xlim,
      "dens" = histogram.coverage,
      "dist" = "coverage"
    )

    # Gene Model
    gene.model = coverage.model.obj$gene.model[[i]]

    # Fitted models & p
    p = point.vals[[i]]
    result = models[[i]]

    # Extract fitted distributions of majority vote model
    majority.vote.plt.data = lapply(result, function(m) {
      mv.model = m[['majority.vote']]
      x = seq(mv.model$seg.start, mv.model$seg.end, by = 1)
      dens = mv.model$dens(x = seq_along(x), mpar = mv.model$par)
      data.frame("x" = x, "dens" = dens, "dist" = mv.model$dist)
    })

    # Coverage plot
    coverage.plt = coverage.plot(
      histogram.coverage = histogram.coverage,
      dist.data = majority.vote.plt.data,
      p = p
    )

    # Residual plot
    residual.plt = residual.plot(
      histogram.coverage = histogram.coverage,
      dist.data = majority.vote.plt.data
    )

    # Metrics plot
    metrics.plt = metrics.plot(
      models = result,
      x.limit = xlim
    )
    min.at = metrics.plt[[3]]
    metrics.plt = metrics.plt[1:2]

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
        colours = c("darkorange", "chartreuse4", "darkorchid4"),
        labels = c("norm", "gamma", "unif"),
        title = expression(bold(underline('Distributions'))),
        lwd = 0.5
      ),
      legend = list(
        colours = c('dodgerblue4', 'gold'),
        labels = c(formatC(min.at, format = "g", digits = 2), 1),
        title = bquote(bold(underline(.("Jaccard Index")))),
        continuous = TRUE,
        height = 2,
        angle = -90,
        tck = 1,
        tck.number = 3,
        at = c(0,100),
        labels.rot = 0,
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

    # Compile plot
    distplot.height = min(2.8, length(result[[1]])*0.8)

    compiled.plts[[i]] = BoutrosLab.plotting.general::create.multipanelplot(
      plot.objects = c(list(coverage.plt, residual.plt), metrics.plt),
      plot.objects.heights = c(10, 3, 1, distplot.height),
      y.spacing = -5,
      # Labels
      main = i,
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

  }

  ## Generate file
  pdf(file.name, width = 10, height = 10)
  print(compiled.plts)
  dev.off()
}
