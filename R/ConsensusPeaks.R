
#' ConsensusPeaks: A wrapper for peak merging functions
#'
#' @param genes A character vector of genes to be tested
#' @param peaks A data frame containing the following columns, and potentially extras, usually found in a BED12 file, base 0 system
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
#' @param rna.or.dna One of 'rna' or 'dna' indicating whether peaks are identified on the transcriptome or genome
#' @param method One 'sf' or 'union'
#' \describe{
#'  \item{sf}{Segment and fit, segmentation of the union set of peaks based on local minima, fitting of continuous distributions}
#'  \item{union}{The union of overlapping peaks is reported}
#' }
#' @param gtf Only if the rna.or.dna parameter is set to 'rna'. The path to a GTF file, character string.
#' @param annotation Only if the rna.or.dna parameter is set to 'rna'. A user can choose to provide a custom annotation file created using the read.gtf function for software efficiency.
#' @param diagnostic Only if the method parameter is set to 'sf'. A logical value indicating whether diagnostic plots for fitted distributions should be plotted.
#' @param fit.mixtures Only if the method parameter is set to 'sf'. A vector of distribution names to be fitted including "unif", "tnorm", "tgamma", "tgamma_flipped", "mixEM" (misxture of normals). "all" for all available distributions.
#' @param trim.step.size Only if the method parameter is set to 'sf'. An integer value > 0 indicating the number of base pairs to trim per iteration from either end of peak to maximize distribution fit quality.
#' @param trim.peak.threshold Only if the method parameter is set to 'sf'. A numeric value between 0 and 1 indicating the maximum proportion of the trimmed peak permitted for distribution fit optimization.
#' @param plot.merged.peaks Only if the method parameter is set to 'sf'. Either a logical value (TRUE or FALSE) indicating all or none of the merged peaks should be plotted. Otherwise, a character vector of genes whose merged peaks should be plotted.
#' @param output.tag A character string added to the names of any output files.
#' @param output.dir Output directory. If the directory does not exist, ConsensusPeaks will attempt to create the directory.
#' @param write.output A logical value indicating whether the output table should be written.
#'
#' @return A data frame with the same columns as 'peaks' indicating the coordinates of the merged peaks in genomic coordinates. If the method is 'dpc', additional columns with the sample names of the samples will indicate the merged p-value using Fisher's Method.
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
#' genes = c("ENSGXX", "ENSGYY"),
#' peaks = peaks,
#' rna.or.dna = "rna",
#' method = "sf",
#' gtf = gtf,
#' annotation = NULL,
#' diagnostic = F,
#' fit.mixture = T,
#' trim.step.size = 10,
#' trim.peak.threshold = 0.1,
#' plot.merged.peaks = F,
#' output.tag = "TEST",
#' output.dir = ".",
#' write.output = F
#' )
#'
ConsensusPeaks = function(
  genes = "all",
  peaks,
  rna.or.dna = c("rna", "dna"),
  method = c("union", "sf"),
  gtf = NULL,
  annotation=NULL,
  diagnostic = F,
  fit.mixtures = "all",
  trim.step.size = 10,
  trim.peak.threshold = 0.1,
  plot.merged.peaks = F,
  output.tag = "",
  output.dir = ".",
  write.output = T
) {

  # Check format of peaks, create list of genes to be tested, create list of plots to be generated
  .check.peaks(peaks)
  if(length(genes) == 1 & genes[1] == "all"){genes = sort(unique(peaks$name))}

  # Check which genes to plot
  if(isTRUE(plot.merged.peaks)){
    plot.merged.peaks = sort(unique(peaks$name))
  } else if(is.logical(plot.merged.peaks)){
    plot.merged.peaks = ""
  } else if(!is.character(plot.merged.peaks)){
    stop("Please provide a character vector of gene names for plot.merged.peaks or a logical T (for all genes) or F (for no genes) to be plotted.")
  }

  # Parameter error checking, generic
  if(!method %in% c("union", "sf")){stop("Please select a method out of 'union' or 'sf'")}
  if(!rna.or.dna %in% c("rna", "dna")){stop("Please select if peaks are part of the transcriptome or genome")}
  if(rna.or.dna == "rna" & (is.null(gtf) & is.null(annotation))){stop("Please provide the GTF file or an annotation to call RNA peaks")}

  # Check the parameters required for sf
  if(method == "sf"){
    if(!is.character(fit.mixtures)){stop("Please provide a character vector for fit.mixture")}
    if(fit.mixtures == "all"){
      fit.mixtures = c("unif", "tnorm", "tgamma", "tgamma_flipped", "mixEM")
    } else {
      fit.mixtures = intersect(fit.mixtures, c("unif", "tnorm", "tgamma", "tgamma_flipped", "mixEM"))
      if(length(fit.mixtures) == 0){stop("Please provide valid distributions")}
    }
    if(!is.logical(diagnostic)){stop("Please provide a logical for diagnostic")}
    if(!is.numeric(trim.step.size) | trim.step.size %% 1 > 0 | trim.step.size < 0){stop("Please provide a positive integer trim.step.size")}
    if(!is.numeric(trim.peak.threshold) | trim.peak.threshold > 1 | trim.peak.threshold < 0){stop("Please provide a numeric between 0 and 1 for trim.peak.threshold")}
  }

  # Error checking, output parameters
  if(!dir.exists(output.dir)){dir.create(output.dir, recursive = T)}
  if(!is.character(output.tag)){stop("Please provide a character output.tag")}
  if(!is.logical(write.output)){stop("Please provide a logical write.output")}

  # All Samples
  all.samples = sort(unique(peaks$sample))

  # If RNA, generate annotation & change peak coordinates
  # TODO: Check annotation file if provided, update function
  if(is.null(annotation) && rna.or.dna == "rna"){annotation = read.gtf(gtf)}

  # Loop through all genes using the correct method
  output.table = data.frame()
  n.genes = length(genes)
  for(i in 1:n.genes){

    # Progress
    i.gene = genes[i]
    i.pct = i/n.genes*100
    print(sprintf("Processing gene: %s", i.gene))
    print(paste0(signif(i.pct, digits = 3), "%"))

    # Methods
    if(method == 'union'){
      results = union.peaks(
        gene = i.gene,
        annotation = annotation,
        peaks = peaks,
        all.samples = all.samples)
    } else if (method == "sf") {
      results = segment.and.fit(
        gene = i.gene,
        annotation = annotation,
        peaks = peaks,
        all.samples = all.samples,
        output.tag = output.tag,
        output.dir = output.dir,
        plot.merged.peaks = plot.merged.peaks,
        diagnostic = diagnostic,
        fit.mixtures = fit.mixtures)
    }
    output.table = rbind(output.table, results)
  }

  # Writing output
  if(write.output){
    filename = file.path(output.dir, paste0(output.tag, ".MergedPeak.tsv"))
    write.table(
      output.table,
      file = filename,
      sep = "\t",
      col.names = T,
      row.names = F,
      quote = F
    )
  }
  # Return output
  return(output.table)
}
