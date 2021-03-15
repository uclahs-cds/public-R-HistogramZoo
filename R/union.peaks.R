union.peaks = function(
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
  PEAKSGR = .retrieve.peaks.as.granges(PEAKS = PEAKS, GENE = GENE, DF = F)
  
  # Get Gene Information
  GENEINFO = .get.gene.anno(GENE, ANNOTATION)
  
  # Converting to RNA
  GENEPEAKSGR = GenomicRanges::shift(PEAKSGR, -1*GENEINFO$left+1)
  GenomicRanges::start(GENEPEAKSGR) = GENEINFO$DNA2RNA[GenomicRanges::start(GENEPEAKSGR)+1]
  GenomicRanges::end(GENEPEAKSGR) = GENEINFO$DNA2RNA[GenomicRanges::end(GENEPEAKSGR)]
  
  # Reduce Overlapping Peaks in the Same Sample
  GENEPEAKSGR = GenomicRanges::reduce(GENEPEAKSGR)
  S4Vectors::mcols(GENEPEAKSGR)$name = GENEINFO$gene
  S4Vectors::mcols(GENEPEAKSGR)$i = 1:length(GENEPEAKSGR)
  
  # Generating Peaks
  merged.peaks.genome = .rna.peaks.to.genome(GENEPEAKSGR, GENEINFO)

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
