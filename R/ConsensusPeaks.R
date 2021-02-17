
ConsensusPeaks = function(
  GENES = "all",
  PEAKS,
  RNA.OR.DNA = c("rna", "dna"),
  METHOD = c("dpc", "union", "corces"),
  GTF = NULL,
  DP.RESOLUTION = 50,
  DP.ITERATIONS = 1000,
  DP.WEIGHT.THRESHOLD = 0.2,
  DP.N.SD = 1,
  DP.ALPHA.PRIORS = c(1,2),
  DP.SEED = 123,
  PLOT.MERGED.PEAKS=F,
  OUTPUT.TAG="",
  OUTPUTDIR =".",
  WRITE.OUTPUT=T

) {

  # Check format of PEAKS, create list of genes to be tested, create list of plots to be generated
  .check.peaks(PEAKS)
  if(GENES == "all"){genes = sort(unique(PEAKS$name))}

  # Check which genes to plot
  if(isTRUE(PLOT.MERGED.PEAKS)){
    PLOT.MERGED.PEAKS = sort(unique(PEAKS$name))
  } else if(isFALSE(PLOT.MERGED.PEAKS)){
    PLOT.MERGED.PEAKS = ""
  } else if(!is.character(PLOT.MERGED.PEAKS)){
    stop("Please provide a character vector of gene names for PLOT.MERGED.PEAKS or a logical T (for all genes) or F (for no genes) to be plotted.")
  }

  # Error checking, generic
  if(!METHOD %in% c("dpc", "union", "corces")){stop("Please select a method out of 'dpc', 'union' or 'corces'")}
  if(!RNA.OR.DNA %in% c("rna", "dna")){stop("Please select if peaks are part of the transcriptome or genome")}
  if(RNA.OR.DNA == "rna" & is.null(GTF)){stop("Please provide the GTF file used to call RNA peaks")}
  # Error checking, dpc
  if(METHOD == "dpc" & !is.integer(DP.RESOLUTION)){stop("Please provide an integer for DP.RESOLUTION")}
  if(METHOD == "dpc" & !is.integer(DP.ITERATIONS)){stop("Please provide an integer for DP.ITERATIONS")}
  if(METHOD == "dpc" & !is.numeric(DP.WEIGHT.THRESHOLD) | DP.WEIGHT.THRESHOLD > 1 | DP.WEIGHT.THRESHOLD < 0){stop("Please provide an number between 0 and 1 for DP.WEIGHT.THRESHOLD")}
  if(METHOD == "dpc" & !is.numeric(DP.N.SD) | DP.N.SD < 0){stop("Please provide an number greater than 0 for DP.N.SD")}
  if(METHOD == "dpc" & !is.integer(DP.SEED)){stop("Please provide an integer for DP.SEED")}
  # Error checking, output
  if(!dir.exists(OUTPUTDIR)){dir.create(OUTPUTDIR, recursive = T)}
  if(!is.character(OUTPUT.TAG)){stop("Please provide a character OUTPUT.TAG")}
  if(!is.logical(WRITE.OUTPUT)){stop("Please provide a logical WRITE.OUTPUT")}

  # Making a list of parameters to pass back and forth
  PARAMETERS = list()
  # General
  PARAMETERS$RNA.OR.DNA = RNA.OR.DNA
  PARAMETERS$METHOD = METHOD
  PARAMETERS$GTF = GTF
  # DP
  PARAMETERS$DP.RESOLUTION = DP.RESOLUTION
  PARAMETERS$DP.ITERATIONS = DP.ITERATIONS
  PARAMETERS$DP.WEIGHT.THRESHOLD = DP.WEIGHT.THRESHOLD
  PARAMETERS$DP.N.SD = DP.N.SD
  PARAMETERS$DP.ALPHA.PRIORS = DP.ALPHA.PRIORS
  PARAMETERS$DP.SEED = DP.SEED
  # Output
  PARAMETERS$OUTPUT.TAG = OUTPUT.TAG
  PARAMETERS$OUTPUTDIR = OUTPUTDIR
  PARAMETERS$PLOT.MERGED.PEAKS = PLOT.MERGED.PEAKS
  # All Samples
  ALL.SAMPLES = sort(unique(PEAKS$sample))
  PARAMETERS$ALL.SAMPLES = ALL.SAMPLES

  # If RNA, generate annotation & change peak coordinates
  ANNOTATION = NULL
  if(PARAMETERS$RNA.OR.DNA == "rna"){
    ANNOTATION = read.gtf(PARAMETERS)
  }
  PEAKSGR = .retrieve.peaks.as.granges(PEAKS = PEAKS, GENE = PARAMETERS$GENE, DF = F)

  # Loop through all genes using the correct method
  OUTPUT.TABLE = data.frame()
  for(i in GENES){
    if(METHOD == 'dpc'){
      RESULTS = dpc(GENE = i, PARAMETERS = PARAMETERS, ANNOTATION = ANNOTATION, PEAKSGR = PEAKSGR, ALL.SAMPLES = ALL.SAMPLES)
    } else if (METHOD == 'coerces'){
      ## todo
    } else if (METHOD == 'union'){
      ## todo
    }
  OUTPUT.TABLE = rbind(OUTPUT.TABLE, RESULTS)
  }

  # Writing output
  if(PARAMETERS$WRITE.OUTPUT){
    filename = paste0(PARAMETERS$OUTPUTDIR, "/", PARAMETERS$GENE, ".", PARAMETERS$OUTPUT.TAG, ".MergedPeaks.tsv")
    write.table(
      OUTPUT.TABLE,
      file = filename,
      sep = "\t",
      col.names = T,
      row.names = F,
      quote = F
    )
  }
  # Return output
  return(OUTPUT.TABLE)
}
