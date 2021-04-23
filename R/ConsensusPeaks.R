
#' ConsensusPeaks: A wrapper for peak merging functions
#'
#' @param GENES A character vector of genes to be tested
#' @param PEAKS A data frame containing the following columns, and potentially extras, usually found in a BED12 file, base 0 system
#' \describe{
#'   \item{chr}{chromosomes, character, same format as those identified in GTF file}
#'   \item{start}{starting position of the peak, integer. base 0}
#'   \item{end}{end position of the peak, integer, base 0}
#'   \item{name}{gene id, character}
#'   \item{score}{p-value associated with the peak}
#'   \item{strand}{strand of the gene, only +, -, * are accepted}
#'   \item{blockCount}{number of segments in the peak, integer}
#'   \item{blockSizes}{size of segments in the peak, BED12 notation, comma-separated}
#'   \item{blockStarts}{starting positions of segments, BED12 notation, comma-separated}
#'   \item{sample}{sample_id of samples, character}
#' }
#' @param RNA.OR.DNA One of 'rna' or 'dna' indicating whether peaks are identified on the transcriptome or genome
#' @param METHOD One 'dpc', 'union', 'corces'
#' \describe{
#'  \item{dpc}{A GMM is fitted on the peaks using a Dirichlet process prior}
#'  \item{union}{The union of overlapping peaks is reported}
#'  \item{coerces}{A 'daisy chain' model described in Coerces et al., 2012, Science, where the most statistically significant peak is chosen}
#' }
#' @param GTF The path to a GTF file, character string.
#' @param DP.RESOLUTION If the METHOD is 'dpc', specify the resolution of the provided PEAKS data.frame. This usually corresponds to the bin width chosen when the peaks are called, a positive integer.
#' @param DP.ITERATIONS If the METHOD is 'dpc', specify the number of iterations used to fit the GMM, a positive integer.
#' @param DP.WEIGHT.THRESHOLD If the METHOD is 'dpc', specify the minimum weight used to threshold the Gaussians. This value is required to be numeric and between 0 and 1.
#' @param DP.N.SD If the METHOD is 'dpc', specify the standard deviation of the Gaussians required to identify peak boundaries, numeric value, above 0.
#' @param DP.ALPHA.PRIORS If the METHOD is 'dpc', specify the priors (alpha, beta) for the Gamma distribution from which the weight concentration prior is drawn, numeric vector of length 2
#' @param DP.SEED A seed for reproducibility, integer above 0.
#' @param PLOT.MERGED.PEAKS Either a logical value (TRUE or FALSE) indicating all or none of the merged peaks should be plotted. Otherwise, a character vector of genes whose merged peaks should be plotted.
#' @param OUTPUT.TAG A character string added to the names of any output files.
#' @param OUTPUTDIR Output directory. If the directory does not exist, ConsensusPeaks will attempt to create the directory.
#' @param WRITE.OUTPUT A logical value indicating whether the output table should be written.
#' @param DIAGNOSTIC
#' @param FIT.MIXTURE
#' @param ANNOTATION
#'
#' @return A data frame with the same columns as 'PEAKS' indicating the coordinates of the merged peaks in genomic coordinates. If the method is 'dpc', additional columns with the sample names of the samples will indicate the merged p-value using Fisher's Method.
#'
#' @export ConsensusPeaks
#'
#' @examples
#' data.path = system.file("extdata", package = "ConsensusPeaks")
#' gtf = paste0(data.path, "/test.gtf")
#' peaks.file = paste0(data.path, "/peaks.bed")
#' peaks = read.table(peaks.file, header = F, sep = '\t', stringsAsFactors = F)
#' colnames(peaks) = c("chr", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRGB", "blockCount", "blockSizes", "blockStarts", "sample")
#'
#' results = ConsensusPeaks(
#' GENES = c("ENSGXX", "ENSGYY"),
#' PEAKS = peaks,
#' RNA.OR.DNA = "rna",
#' METHOD = "dpc",
#' GTF = gtf,
#' DP.RESOLUTION = 10,
#' DP.ITERATIONS = 1000,
#' DP.WEIGHT.THRESHOLD = 0.1,
#' DP.N.SD = 1.5,
#' DP.ALPHA.PRIORS = c(0.5,5),
#' DP.SEED = 123,
#' PLOT.MERGED.PEAKS = F,
#' OUTPUT.TAG = "TEST",
#' OUTPUTDIR = ".",
#' WRITE.OUTPUT = F)
#'
ConsensusPeaks = function(
  GENES = "all",
  PEAKS,
  RNA.OR.DNA = c("rna", "dna"),
  METHOD = c("union", "sf"),
  GTF = NULL,
  PLOT.MERGED.PEAKS=F,
  OUTPUT.TAG="",
  OUTPUTDIR =".",
  WRITE.OUTPUT=T,
  DIAGNOSTIC = F,
  FIT.MIXTURE = T,
  ANNOTATION=NULL
) {

  # Check format of PEAKS, create list of genes to be tested, create list of plots to be generated
  .check.peaks(PEAKS)
  if(length(GENES) == 1 & GENES[1] == "all"){genes = sort(unique(PEAKS$name))}

  # Check which genes to plot
  if(isTRUE(PLOT.MERGED.PEAKS)){
    PLOT.MERGED.PEAKS = sort(unique(PEAKS$name))
  } else if(is.logical(PLOT.MERGED.PEAKS)){
    PLOT.MERGED.PEAKS = ""
  } else if(!is.character(PLOT.MERGED.PEAKS)){
    stop("Please provide a character vector of gene names for PLOT.MERGED.PEAKS or a logical T (for all genes) or F (for no genes) to be plotted.")
  }

  # Error checking, generic
  if(!METHOD %in% c("union", "sf")){stop("Please select a method out of 'union' or 'sf'")}
  if(!RNA.OR.DNA %in% c("rna", "dna")){stop("Please select if peaks are part of the transcriptome or genome")}
  if(RNA.OR.DNA == "rna" & is.null(GTF)){stop("Please provide the GTF file used to call RNA peaks")}
  # Error checking, output
  if(!dir.exists(OUTPUTDIR)){dir.create(OUTPUTDIR, recursive = T)}
  if(!is.character(OUTPUT.TAG)){stop("Please provide a character OUTPUT.TAG")}
  if(!is.logical(WRITE.OUTPUT)){stop("Please provide a logical WRITE.OUTPUT")}

  # All Samples
  ALL.SAMPLES = sort(unique(PEAKS$sample))

  # If RNA, generate annotation & change peak coordinates
  if(is.null(ANNOTATION) && RNA.OR.DNA == "rna"){
    ANNOTATION = read.gtf(GTF)
  }

  # Loop through all genes using the correct method
  OUTPUT.TABLE = data.frame()
  for(i in GENES){
    print(sprintf("Processing gene: %s", i))
    print( paste(as.character(signif(which(GENES == i)/length(GENES)*100, digits = 3)),"%") )
    if(METHOD == 'union'){
      RESULTS = union.peaks(GENE = i, ANNOTATION = ANNOTATION, PEAKS = PEAKS, ALL.SAMPLES = ALL.SAMPLES)
    } else if (METHOD == "sf"){
      RESULTS = segment.and.fit(
        GENE = i,
        ANNOTATION = ANNOTATION,
        PEAKS = PEAKS,
        ALL.SAMPLES = ALL.SAMPLES,
        OUTPUT.TAG = OUTPUT.TAG,
        OUTPUTDIR = OUTPUTDIR,
        PLOT.MERGED.PEAKS = PLOT.MERGED.PEAKS,
        DIAGNOSTIC = DIAGNOSTIC,
        FIT.MIXTURE = FIT.MIXTURE)
    }
  OUTPUT.TABLE = rbind(OUTPUT.TABLE, RESULTS)
  }

  # Writing output
  if(WRITE.OUTPUT){
    filename = file.path(OUTPUTDIR, paste0(OUTPUT.TAG, ".MergedPeak.tsv"))
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
