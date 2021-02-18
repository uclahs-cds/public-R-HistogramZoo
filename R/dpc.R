dpc = function(
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

  # Examining Weights & Distributions of the GMM
  dp = ConsensusPeaks:::.dpc.peaks(GENEPEAKSGR, PARAMETERS)

  # Generating Peaks
  merged.peaks.rna = ConsensusPeaks:::.generate.peaks.from.gmm(dp = dp, PARAMETERS, GENEINFO)
  merged.peaks.genome = ConsensusPeaks:::.rna.peaks.to.genome(merged.peaks.rna, GENEINFO)

  # Plotting
  if(GENE %in% PARAMETERS$PLOT.MERGED.PEAKS){
    plotting.data = ConsensusPeaks:::.generate.merged.peaks.plotting(dp, PARAMETERS, GENEINFO, GENEPEAKSGR)
    ConsensusPeaks:::.plot.merged.peaks(GENE, plotting.data, merged.peaks.rna, PARAMETERS)
  }

  # Return a Data Frame of Merged Peaks
  if(length(merged.peaks.genome) == 0){
    warning("No Peaks are Found for This Gene After Fitting!", call. = TRUE, domain = NULL)
    OUTPUT.TABLE = .generate.null.result(PARAMETERS)
  } else {

    # Creating a BED12 File
    GenomicRanges::start(merged.peaks.genome) = GenomicRanges::start(merged.peaks.genome)-1
    PEAKS.FINAL = .bed6tobed12(MERGED.PEAKS = merged.peaks.genome, ID.COLS = c("name", "i"))

    # Merging P-Values
    SAMPLE.PVAL = .merge.p(PEAKSGR, MERGED.PEAKS = merged.peaks.genome, ANNOTATION, PARAMETERS, ID.COLS = c("name", "i"))

    # Write Output Tables & Return Files
    OUTPUT.TABLE = merge(PEAKS.FINAL, SAMPLE.PVAL, by = "peak", all = T)
    OUTPUT.TABLE = OUTPUT.TABLE[,colnames(OUTPUT.TABLE) != "peak"]
  }

  return(OUTPUT.TABLE)

}
