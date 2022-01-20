
coverage.plot = function(
  histogram.coverage,
  p = NULL
){
  
  xlim = length(histogram.coverage)
  plt.data = data.frame(
    "x" = 1:xlim,
    "dens" = histogram.coverage,
    "dist" = "coverage",
    "i" = 1
  )
  points = plt.data[plt.data$x %in% p,]
  
  # Reference Colours, Numeric & Line Width for Distributions
  col.reference = c("black", "darkorange", "chartreuse4", "darkorchid4")
  names(col.reference) = c("coverage", "norm", "gamma", "unif")
  lwd.reference = c(1, 2.5, 2.5, 2.5)
  names(lwd.reference) = c("coverage", "norm", "gamma", "unif")
  
  # Plotting attributes for lineplot
  col.used = unique(lineplot.data[,c("dist", "i")])
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

plot.fitted.segments = function(
  histogram.names,
  coverage.model.obj
){
 
  # Loop through genes to plot things
  for(i in histogram.names){
    
    histogram.coverage = coverage.model.obj$histogram.coverage[[i]]
    gene.model = coverage.model.obj$gene.model[[i]]
    result = coverage.model.obj$results[[i]]
    
    # Coverage plot
    # We need p and a way of generating the distribution
    coverage.plt = coverage.plot(
      histogram.coverage = histogram.coverage,
      p = NULL
    )
    
    # Residual plot
    
    # Peaks and metrics
    
    # Transcript plot
    if(!is.null(gene.model)){
      
    }
    
    # Legend
    # Compile plot
    
    
  }
  
  ## Generate file 
}