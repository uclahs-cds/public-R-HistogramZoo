fit.residuals.threshold.unif = function(
  x,
  mod,
  plot.diagnostic.residuals = F,
  output.dir,
  i,
  gene
){

  # Initializing variables
  bin.data = table(x)
  bin.count = as.integer(bin.data)
  bin.x = as.numeric(names(bin.data))
  scalefactor = length(x)

  # Observed
  density.table = data.frame(x = bin.x, Observed = bin.count)

  # Expected
  dens.scale = calculate.density(
    m = m,
    x = bin.x,
    seg.start = NULL,
    seg.end = NULL,
    stepsize = 1,
    scale.density = T,
    return.df = F)

  # Residuals
  fit.residuals = dens.scale - bin.count
  fit.residuals = fit.residuals/scalefactor

  # Plotting Data
  density.table[,'Expected'] = dens.scale
  density.table[,'Residuals'] = fit.residuals


  if(plot.diagnostic.residuals){

    filename = file.path(output.dir, paste0(gene, ".residuals.uniform.", i, ".pdf"))
    pdf(filename)

    p1 = ggplot2::ggplot(density.table, ggplot2::aes(x = x, y = Observed)) +
      ggplot2::geom_step() +
      ggplot2::geom_step(ggplot2::aes(y = Expected, col = 'darkorchid4')) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::ggtitle(paste0(gene, ": ", i)) +
      ggplot2::ylab("Coverage (at BP resolution)") +
      ggplot2::xlab("Scaled X")

    p2 = ggplot2::ggplot(density.table, ggplot2::aes(x = x, y = Residuals)) +
      ggplot2::geom_point() +
      ggplot2::theme_bw() +
      ggplot2::ggtitle("Uniform Fit Residuals")

    # Printing the plot
    print(p1)
    print(p2)

    dev.off()
  }

  return(density.table)

}
