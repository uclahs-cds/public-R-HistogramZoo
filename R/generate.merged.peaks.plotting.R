.generate.merged.peaks.plotting = function(dp, PARAMETERS, GENEINFO, GENEPEAKSGR){

  # Peak Coverage
  REDUCED.GENE.PEAKS.GR = GenomicRanges::reduce(GENEPEAKSGR)
  PEAK.COVERAGE = GenomicRanges::coverage(GENEPEAKSGR)
  BINS = unlist(GenomicRanges::tile(REDUCED.GENE.PEAKS.GR, width = 1))
  BIN.COUNTS = data.frame(GenomicRanges::binnedAverage(BINS, PEAK.COVERAGE, "Coverage"), stringsAsFactors = F)

  # Plotting Data
  x.norm <- seq(min(dp$startvec.scaled), max(dp$startvec.scaled), by=0.01)
  y.fit <- data.frame(replicate(100, dirichletprocess::PosteriorFunction(dp$dp)(x.norm)))
  fit.frame <- data.frame(x=x.norm, y=rowMeans(y.fit))
  fit.frame$x = (fit.frame$x*dp$startvec.sd)+dp$startvec.mean
  fit.frame$y = fit.frame$y*max(BIN.COUNTS$Coverage)

  # Plotting DP data
  dp_data = dp[['dp_data']]
  sample.points = seq(1, GENEINFO$exome_length, 10)
  scaling.factor = length(PARAMETERS$ALL.SAMPLES)*100
  plot.dp.data = data.frame("sample.points" = sample.points, stringsAsFactors = F)
  for(i in 1:nrow(dp_data)){
    norm.tmp = dnorm(sample.points, mean = dp_data$Mu[i], sd = dp_data$Sigma[i])*dp_data$Weights[i]*scaling.factor
    plot.dp.data = cbind(plot.dp.data, norm.tmp)
  }
  colnames(plot.dp.data) = c("sample.points", paste0("V", 1:nrow(dp_data)))

  list("bin.counts" = BIN.COUNTS, "fit.frame" = fit.frame, "plot.dp.data" = plot.dp.data, "startvec" = dp$startvec)
}
