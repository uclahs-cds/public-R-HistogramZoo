
ggplot.plot = function(
  outputdir,
  distr.plotting.data,
  geneinfo,
  bin.counts,
  seg.gr,
  p
) {

  filename = file.path(outputdir, paste0(geneinfo$gene, ".SegmentAndFit.pdf"))
  pdf(filename, width = 10, height = 10)

  p1 = ggplot2::ggplot(bin.counts, ggplot2::aes(x = start, y = Coverage)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(data = bin.counts[p,], ggplot2::aes(x = start, y = Coverage), col = 'red') +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(geneinfo$gene) +
    ggplot2::ylab("Coverage (at BP resolution)") +
    ggplot2::xlab("Transcript Coordinate") +
    ggplot2::annotate("rect", xmin=GenomicRanges::start(seg.gr), xmax=GenomicRanges::end(seg.gr), ymin=-1 , ymax=-0.1, alpha=0.5, color="black", fill=1:length(seg.gr))

  for(i in 1:length(distr.plotting.data)){
    p1 = p1 + ggplot2::geom_line(data = distr.plotting.data[[i]], ggplot2::aes(x=x, y=dens, color = col))
  }
  p1 = p1 + ggplot2::guides(col=ggplot2::guide_legend(title="Distribution"))
  print(p1)

  dev.off()
}
