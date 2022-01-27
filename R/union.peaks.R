union.peaks = function(
  gene,
  annotation,
  peaks,
  all.samples
){

  # If the gene doesn't have peaks
  if(!gene %in% peaks$name){
    warn.message = paste0("No Peaks are Found for ", gene, " in peaks!")
    warning(warn.message, call. = TRUE, domain = NULL)
    return(.generate.null.result(all.samples))
  }

  # peaksgr
  peaksgr = .retrieve.peaks.as.granges(peaks = peaks, gene = gene, return.df = F)

  # Get Gene Information
  geneinfo = .get.gene.anno(gene, annotation)

  # Converting to RNA
  genepeaksgr = GenomicRanges::shift(peaksgr, -1*geneinfo$left+1)
  GenomicRanges::start(genepeaksgr) = geneinfo$DNA2RNA[GenomicRanges::start(genepeaksgr)+1]
  GenomicRanges::end(genepeaksgr) = geneinfo$DNA2RNA[GenomicRanges::end(genepeaksgr)]

  # Reduce Overlapping Peaks in the Same Sample
  genepeaksgr = GenomicRanges::reduce(genepeaksgr)
  S4Vectors::mcols(genepeaksgr)$name = geneinfo$gene
  S4Vectors::mcols(genepeaksgr)$i = 1:length(genepeaksgr)

  # Generating Peaks
  merged.peaks.genome = .rna.peaks.to.genome(genepeaksgr, geneinfo)

  # Return a Data Frame of Merged Peaks
  if(length(merged.peaks.genome) == 0){
    warning("No Peaks are Found for This Gene After Fitting!", call. = TRUE, domain = NULL)
    output.table = .generate.null.result(all.samples)
  } else {

    # Creating a BED12 File
    GenomicRanges::start(merged.peaks.genome) = GenomicRanges::start(merged.peaks.genome)-1
    peaks.final = .bed6tobed12(merged.peaks = merged.peaks.genome, id.cols = c("name", "i"))

    # Merging P-Values
    sample.pval = .merge.p(peaksgr, merged.peaks = merged.peaks.genome, annotation, all.samples)

    # Write Output Tables & Return Files
    output.table = merge(peaks.final, sample.pval, by = "peak", all = T)
    output.table = output.table[,colnames(output.table) != "peak"]
  }

  return(output.table)

}
