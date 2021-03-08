
hmm = function(
  GENE,
  PARAMETERS,
  ANNOTATION,
  PEAKS
){

  # PEAKSGR
  PEAKSGR = .retrieve.peaks.as.granges(PEAKS = PEAKS, GENE = GENE, DF = F)

  # Get Gene Information
  GENEINFO = .get.gene.anno(GENE, ANNOTATION)

  # Converting to RNA
  GENEPEAKSGR = GenomicRanges::shift(PEAKSGR, -1*GENEINFO$left+1)
  GenomicRanges::start(GENEPEAKSGR) = GENEINFO$DNA2RNA[GenomicRanges::start(GENEPEAKSGR)+1]
  GenomicRanges::end(GENEPEAKSGR) = GENEINFO$DNA2RNA[GenomicRanges::end(GENEPEAKSGR)]

  # Reduce Overlapping Peaks in the Same Sample
  GENEPEAKSGR = S4Vectors::split(GENEPEAKSGR, GENEPEAKSGR$sample)
  GENEPEAKSGR = unlist(GenomicRanges::reduce(GENEPEAKSGR))

  # Tiling Peaks
  TILED.PEAKS.GR = unlist(GenomicRanges::tile(GENEPEAKSGR, width = PARAMETERS$DP.RESOLUTION))

  # Peak Coverage
  REDUCED.GENE.PEAKS.GR = GenomicRanges::reduce(GENEPEAKSGR)
  PEAK.COVERAGE = GenomicRanges::coverage(GENEPEAKSGR)
  BINS = unlist(GenomicRanges::tile(REDUCED.GENE.PEAKS.GR, width = 1))
  BIN.COUNTS = data.frame(GenomicRanges::binnedAverage(BINS, PEAK.COVERAGE, "Coverage"), stringsAsFactors = F)
  add.zero.counts = setdiff(1:GENEINFO$exome_length, BIN.COUNTS$start)
  add.zero.counts = data.frame("start" = add.zero.counts, "Coverage" = 0, stringsAsFactors = F)
  BIN.COUNTS = rbind(BIN.COUNTS[,c("start", "Coverage")], add.zero.counts)
  BIN.COUNTS = BIN.COUNTS[order(BIN.COUNTS$start),]


  # HMM
  result = tryCatch({
    # HMM
    mod <- depmixS4::depmix(Coverage ~ 1, data = BIN.COUNTS, nstates = 2, family = gaussian())
    fit.mod <- depmixS4::fit(mod)
    est.states <- depmixS4::posterior(fit.mod)
    peak.state = ifelse(mean(BIN.COUNTS$Coverage[est.states$state == 1]) > mean(BIN.COUNTS$Coverage[est.states$state == 2]), 1, 2)

    # Generating Peaks
    rna.peaks.gr = GenomicRanges::GRanges(seqnames = GENEINFO$chr,
                                          IRanges::IRanges(which(est.states$state == peak.state)),
                                          strand = GENEINFO$strand)

  }, finally = {
    # Taking the Union of Peaks Regions if HMM fails to fit
    rna.peaks.gr = GenomicRanges::GRanges(seqnames = GENEINFO$chr,
                                          IRanges::IRanges(which(BIN.COUNTS$Coverage > 0)),
                                          strand = GENEINFO$strand)
  })
  mcols(rna.peaks.gr)$name = GENEINFO$gene
  mcols(rna.peaks.gr)$i = 1:length(rna.peaks.gr)

  # DNA peak
  dna.peaks.gr = .rna.peaks.to.genome(rna.peaks.gr, GENEINFO)

  # Plotting
  if(GENE %in% PARAMETERS$PLOT.MERGED.PEAKS){

    filename = paste0(PARAMETERS$OUTPUTDIR, "/", GENE, ".", PARAMETERS$OUTPUT.TAG, ".MergedPeaks.pdf")
    pdf(filename)

    p2 = ggplot2::ggplot() +
      ggplot2::geom_line(data=BIN.COUNTS, ggplot2::aes(x=start, y=Coverage), colour='grey') +
      ggplot2::theme_bw() +
      ggplot2::ggtitle(GENE) +
      ggplot2::ylab("Coverage (at BP resolution)") +
      ggplot2::xlab("Transcript Coordinate") +
      ggplot2::annotate("rect", xmin=GenomicRanges::start(rna.peaks.gr), xmax=GenomicRanges::end(rna.peaks.gr), ymin=-1 , ymax=-0.1, alpha=0.5, color="black", fill=length(rna.peaks.gr))

    print(p2)
    dev.off()
  }

  # Return a Data Frame of Merged Peaks
  if(length(dna.peaks.gr) == 0){
    warning("No Peaks are Found for This Gene After Fitting!", call. = TRUE, domain = NULL)
    OUTPUT.TABLE = .generate.null.result(PARAMETERS)
  } else {

    # Creating a BED12 File
    GenomicRanges::start(dna.peaks.gr) = GenomicRanges::start(dna.peaks.gr)-1
    PEAKS.FINAL = .bed6tobed12(MERGED.PEAKS = dna.peaks.gr, ID.COLS = c("name", "i"))

    # Merging P-Values
    SAMPLE.PVAL = .merge.p(PEAKSGR, MERGED.PEAKS = dna.peaks.gr, ANNOTATION, PARAMETERS, ID.COLS = c("name", "i"))

    # Write Output Tables & Return Files
    OUTPUT.TABLE = merge(PEAKS.FINAL, SAMPLE.PVAL, by = "peak", all = T)
    OUTPUT.TABLE = OUTPUT.TABLE[,colnames(OUTPUT.TABLE) != "peak"]
  }

  return(OUTPUT.TABLE)

}
