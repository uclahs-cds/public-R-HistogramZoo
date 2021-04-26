fit.residuals = function(
  x,
  mod,
  plot.diagnostic.residuals = F,
  output.dir,
  i,
  gene
){

  # Observed
  bin.data = table(x)
  scalefactor = length(x)

  # Residual Table
  bin.count = as.integer(bin.data)
  # x value
  bin.x = as.numeric(names(bin.data))

  density.table = data.frame(x = bin.x, Observed = bin.count)

  abs.residual.table = data.frame("x" = bin.x)
  sq.residual.table = abs.residual.table
  residual.table = abs.residual.table
  density.list = list()

  # Colour
  col.reference = structure(c("black", "orange", "chartreuse4", "chartreuse3", "darkorange", "darkorchid4"),
                            names = c("coverage", "tnorm", "tgamma", "tgamma_flip", "norm_mixture", "unif"))

  # Expected
  for(m in mod){

    dens.scale = calculate.density(
      m = m,
      x = bin.x,
      seg.start = NULL,
      seg.end = NULL,
      stepsize = 1,
      scale.density = T,
      return.df = F)

    fit.residuals = dens.scale - bin.count
    fit.residuals = fit.residuals/scalefactor
    # fit.residuals <- (dens.scale - bin.count)^2

    # Adding to table
    distname = ifelse(class(m) == "mixEM", "mixEM", m$distname)
    residual.table[,distname] = fit.residuals
    sq.residual.table[,distname] = fit.residuals^2
    abs.residual.table[,distname] = abs(fit.residuals)
    density.list[[distname]] = data.frame("x" = bin.x, "y" = dens.scale, "col" = distname, stringsAsFactors = F)
  }

  if(plot.diagnostic.residuals){

    filename = file.path(output.dir, paste0(gene, ".residuals.", i, ".pdf"))
    pdf(filename)

    plotting.data = do.call(rbind.data.frame, density.list)
    plotting.data.res = reshape2::melt(residual.table, id.vars = "x")
    plotting.data.sq.res = reshape2::melt(sq.residual.table, id.vars = "x")
    plotting.data.abs.res = reshape2::melt(abs.residual.table, id.vars = "x")
    plotting.data.mean.res = rowMeans(sq.residual.table[,2:ncol(sq.residual.table)])
    plotting.data.mean.res = data.frame("x" = sq.residual.table$x, "MeanModelResidual" = plotting.data.mean.res)

    # Keep same factor level
    plotting.data$col <- factor(plotting.data$col, levels(as.factor(plotting.data.res$variable)))

    p1 = ggplot2::ggplot(density.table, ggplot2::aes(x = x, y = Observed)) +
      ggplot2::geom_step() +
      ggplot2::theme_bw() +
      ggplot2::ggtitle(paste0(gene, "\tSegment:", i)) +
      ggplot2::ylab("Coverage (at BP resolution)") +
      ggplot2::xlab("Scaled X")

    p2 = ggplot2::ggplot(plotting.data.res, ggplot2::aes(x = x, y = value, color = variable)) +
      ggplot2::geom_point() +
      ggplot2::theme_bw() +
      ggplot2::ggtitle("Residuals")

    p3 = ggplot2::ggplot(plotting.data.sq.res, ggplot2::aes(x = x, y = value, color = variable)) +
      ggplot2::geom_point() +
      ggplot2::theme_bw() +
      ggplot2::ggtitle("Squared Residuals")

    p4 = ggplot2::ggplot(plotting.data.abs.res, ggplot2::aes(x = x, y = value, color = variable)) +
      ggplot2::geom_point() +
      ggplot2::theme_bw() +
      ggplot2::ggtitle("Absolute Residuals")

    # Mean Squared Residuals
    p5 = ggplot2::ggplot(plotting.data.mean.res, ggplot2::aes(x = x, y = MeanModelResidual)) +
      ggplot2::geom_line() +
      ggplot2::theme_bw() +
      ggplot2::ggtitle("Mean Model Residual")

    # Printing the plot
    print(p1)
    print(p2)
    print(p3)
    print(p4)
    print(p5)

    dev.off()
  }


}
