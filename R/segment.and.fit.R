
segment.and.fit = function(
  GENE,
  PARAMETERS,
  ANNOTATION,
  PEAKS
){
  
  # If the gene doesn't have peaks
  if(!GENE %in% PEAKS$name){
    warn.message = paste0("No Peaks are Found for ", GENE, " in PEAKS!")
    warning(warn.message, call. = TRUE, domain = NULL)
    return(.generate.null.result(PARAMETERS))
  }
  
  # PEAKSGR
  PEAKSGR = ConsensusPeaks:::.retrieve.peaks.as.granges(PEAKS = PEAKS, GENE = GENE, DF = F)
  
  # Get Gene Information
  GENEINFO = ConsensusPeaks:::.get.gene.anno(GENE, ANNOTATION)
  
  # Converting to RNA
  GENEPEAKSGR = GenomicRanges::shift(PEAKSGR, -1*GENEINFO$left+1)
  GenomicRanges::start(GENEPEAKSGR) = GENEINFO$DNA2RNA[GenomicRanges::start(GENEPEAKSGR)+1]
  GenomicRanges::end(GENEPEAKSGR) = GENEINFO$DNA2RNA[GenomicRanges::end(GENEPEAKSGR)]
  
  # Reduce Overlapping Peaks in the Same Sample
  GENEPEAKSGR = S4Vectors::split(GENEPEAKSGR, GENEPEAKSGR$sample)
  GENEPEAKSGR = unlist(GenomicRanges::reduce(GENEPEAKSGR))
  
  # Peak Coverage
  REDUCED.GENE.PEAKS.GR = GenomicRanges::reduce(GENEPEAKSGR)
  PEAK.COVERAGE = GenomicRanges::coverage(GENEPEAKSGR)
  BINS = unlist(GenomicRanges::tile(REDUCED.GENE.PEAKS.GR, width = 1))
  BIN.COUNTS = data.frame(GenomicRanges::binnedAverage(BINS, PEAK.COVERAGE, "Coverage"), stringsAsFactors = F)
  add.zero.counts = setdiff(1:GENEINFO$exome_length, BIN.COUNTS$start)
  add.zero.counts = data.frame("start" = add.zero.counts, "Coverage" = 0, stringsAsFactors = F)
  BIN.COUNTS = rbind(BIN.COUNTS[,c("start", "Coverage")], add.zero.counts)
  BIN.COUNTS = BIN.COUNTS[order(BIN.COUNTS$start),]
  
  # Segmenting & Determining which segments are peaks
  p <- find.peaks(-BIN.COUNTS$Coverage, m = 150)
  p = sort(p)
  seg.df = data.frame("start" = p[1:(length(p)-1)]+1, "end" = p[2:length(p)], "mean" = 0)
  for(i in 1:nrow(seg.df)){
    tmp = BIN.COUNTS$Coverage[BIN.COUNTS$start >= seg.df$start[i] & BIN.COUNTS$start <= seg.df$end[i]]
    seg.df$mean[i] = mean(tmp)
  }
  seg.df = seg.df[seg.df$mean > 0,]
  
  # Tiling Peaks
  peak.counts = unlist(GenomicRanges::tile(GENEPEAKSGR, width = 1))
  peak.counts = GenomicRanges::start(peak.counts)
  
  # Segments
  filename = "~/figures/segments.pdf"
  pdf(filename, width = 5, height = 5)
  plot(BIN.COUNTS$start, BIN.COUNTS$Coverage, type = "s")
  p <- find.peaks(-BIN.COUNTS$Coverage, m = 150)
  points(BIN.COUNTS$start[p], BIN.COUNTS$Coverage[p], col = 'red')
  dev.off()
  
  # Fitting different models
  results = data.frame()
  for(i in 1:nrow(seg.df)){
    x = peak.counts[peak.counts >= seg.df$start[i] & peak.counts <= seg.df$end[i]]
    
    models = lapply(c("norm", "gamma", "unif"), function(distr) fitdistrplus::fitdist(x, distr))
    fits = lapply(models, function(mod){
      # c("dist" = summary(mod)$distname, "aic" = summary(mod)$aic, "bic" = summary(mod)$bic, summary(mod)$estimate)
      data.frame(
        "dist" = summary(mod)$distname,
        "loglikelihood" = summary(mod)$loglik,
        "aic" = summary(mod)$aic,
        "bic" = summary(mod)$bic
      )
    })
    fits = do.call(rbind, fits)
    fits$i = i
    results = rbind(results, fits)
    
    filename = paste0("~/figures/fit.segments.", i, ".pdf")
    pdf(filename, width = 10, height = 10)
    par(mfrow = c(2,2))
    plot.legend = c("norm", "gamma", "uniform")
    denscomp(models, legendtext = plot.legend)
    qqcomp(models, legendtext = plot.legend)
    cdfcomp(models, legendtext = plot.legend)
    ppcomp(models, legendtext = plot.legend)
    dev.off()
    
  }
  write.table(
    results,
    file = "~/figures/results.fit.segments.tsv",
    sep = "\t",
    col.names = T,
    row.names = F,
    quote = F
  )
  
  
  
}