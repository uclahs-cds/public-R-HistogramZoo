fit.residuals = function(
  x,
  scalefactor,
  mod,
  sample.size,
  plot.diagnostic.residuals = F,
  outputdir,
  i,
  gene
){

  # Observed
  bin.data <- table(x)

  # Residual Table
  # "Observed" = as.vector(bin.data)
  # "x" = as.numeric(names(bin.data)),
  abs.residual.table = data.frame("x" = as.numeric(names(bin.data)))
  sq.residual.table = abs.residual.table
  residual.table = abs.residual.table
  density.list = list()

  # Colour
  col.reference = structure(c("black", "orange", "chartreuse4", "chartreuse3", "darkorange", "darkorchid4"),
                            names = c("coverage", "tnorm", "tgamma", "tgamma_flip", "norm_mixture", "unif"))

  # Expected
  for(m in mod){

    if(class(m) == "mixEM") {
      distname <- "norm_mixture"
      dens <- dnorm_mixture(as.numeric(names(bin.data)), m)
    } else {
      params <- c(as.list(m$estimate), as.list(m$fix.arg))
      distname <- m$distname
      ddistname <- paste0("d", distname)
      call.params <- c(list(x = as.numeric(names(bin.data))), as.list(params))
      dens <- do.call(ddistname, call.params)
    }
    dens.scale <- dens * scalefactor
    fit.residuals = dens.scale - as.integer(bin.data)
    fit.residuals = fit.residuals/sample.size
    # fit.residuals <- (dens.scale - as.integer(bin.data))^2

    # Adding to table
    residual.table[,distname] = fit.residuals
    sq.residual.table[,distname] = fit.residuals^2
    abs.residual.table[,distname] = abs(fit.residuals)
    density.list[[distname]] = data.frame("x" = as.numeric(names(bin.data)), "y" = dens.scale, "col" = as.vector(col.reference[distname]), stringsAsFactors = F)
  }

  if(plot.diagnostic.residuals){

    filename = file.path(outputdir, paste0(GENE, ".residuals.", i, ".pdf"))
    pdf(filename)

    p1 = ggplot2::ggplot(density.table, ggplot2::aes(x = x, y = Observed)) +
      ggplot2::geom_line() +
      ggplot2::theme_bw() +
      ggplot2::ggtitle(gene) +
      ggplot2::ylab("Coverage (at BP resolution)") +
      ggplot2::xlab("Scaled X")

    for(j in 1:length(density.list)){
      p1 = p1 + ggplot2::geom_line(data = density.list[[j]], ggplot2::aes(x=x, y=y, color = col))
    }
    p1 = p1 + ggplot2::guides(col=ggplot2::guide_legend(title="Distribution"))

    plotting.data = reshape2::melt(residual.table, id.vars = "x")
    p2 = ggplot2::ggplot(plotting.data, aes(x = x, y = value, color = variable)) +
      ggplot2::geom_point() +
      ggplot2::theme_bw() +
      ggplot2::ggtitle("Residuals")

    plotting.data = reshape2::melt(sq.residual.table, id.vars = "x")
    p3 = ggplot2::ggplot(plotting.data, aes(x = x, y = value, color = variable)) +
      ggplot2::geom_point() +
      ggplot2::theme_bw() +
      ggplot2::ggtitle("Squared Residuals")

    plotting.data = reshape2::melt(abs.residual.table, id.vars = "x")
    p4 = ggplot2::ggplot(plotting.data, aes(x = x, y = value, color = variable)) +
      ggplot2::geom_point() +
      ggplot2::theme_bw() +
      ggplot2::ggtitle("Absolute Residuals")

   # Mean Squared Residuals
   plotting.data = rowMeans(sq.residual.table[,2:ncol(sq.residual.table)])
   plotting.data = data.frame("x" = sq.residual.table$x, "MeanModelResidual" = plotting.data)
   p5 = ggplot2::ggplot(plotting.data, aes(x = x, y = MeanModelResidual)) +
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
