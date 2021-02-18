.plot.merged.peaks = function(
  GENE,
  plotting.data,
  plot.merged.peaks,
  PARAMETERS
){

  melted.plot.dp.data = reshape2::melt(plotting.data$plot.dp.data, id = "sample.points")

  filename = paste0(PARAMETERS$OUTPUTDIR, "/", GENE, ".", PARAMETERS$OUTPUT.TAG, ".MergedPeaks.pdf")
  pdf(filename)

  p1 = ggplot2::ggplot() +
    ggplot2::geom_histogram(data=data.frame("start" = plotting.data$startvec), ggplot2::aes(x=start), binwidth = 50, colour = "grey", fill="white") +
    ggplot2::geom_line(data=plotting.data$fit.frame, ggplot2::aes(x=x, y=y), colour='black') +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(GENE) +
    ggplot2::ylab("Binned Counts From Peaks (Used to Fit GMM)") +
    ggplot2::xlab("Transcript Coordinate") +
    ggplot2::annotate("rect", xmin=GenomicRanges::start(plot.merged.peaks), xmax=GenomicRanges::end(plot.merged.peaks), ymin=-1 , ymax=-0.1, alpha=0.5, color="black", fill=length(plot.merged.peaks))

  p1 = p1 + ggplot2::geom_line(data = melted.plot.dp.data, ggplot2::aes(x=sample.points, y=value, color = variable)) +
    ggplot2::theme(legend.position = "none")

  print(p1)

  p2 = ggplot2::ggplot() +
    ggplot2::geom_line(data=plotting.data$bin.counts, ggplot2::aes(x=start, y=Coverage), colour='grey') +
    ggplot2::geom_line(data=plotting.data$fit.frame, ggplot2::aes(x=x,y=y), colour='black') +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(GENE) +
    ggplot2::ylab("Coverage (at BP resolution)") +
    ggplot2::xlab("Transcript Coordinate") +
    ggplot2::annotate("rect", xmin=GenomicRanges::start(plot.merged.peaks), xmax=GenomicRanges::end(plot.merged.peaks), ymin=-1 , ymax=-0.1, alpha=0.5, color="black", fill=length(plot.merged.peaks))

  p2 = p2 + ggplot2::geom_line(data = melted.plot.dp.data, ggplot2::aes(x=sample.points, y=value, color = variable)) +
    ggplot2::theme(legend.position = "none")

  print(p2)

  dev.off()

}
