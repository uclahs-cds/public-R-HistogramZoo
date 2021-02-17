
#' Creates a pile up of bed files using the Sushi R package
#'
#' @param GENE A 'character' gene id corresponding to the gene_id's found in the GTF files and the 'name' column in the peak files
#' @param PEAKS data frame containing the following columns, and potentially extras, usually found in a BED12 file, base 0 system
#' \describe{
#'   \item{chr}{chromosomes, same as in GTF file}
#'   \item{start}{starting position of the peak, base 0}
#'   \item{end}{end position of the peak, base 0}
#'   \item{name}{gene id}
#'   \item{score}{p-value associated with the peak}
#'   \item{strand}{strand of the gene}
#'   \item{blockCount}{number of segments in the peak}
#'   \item{blockSizes}{size of segments in the peak, BED12 notation}
#'   \item{blockStarts}{starting positions of segments, BED12 notation}
#'   \item{sample}{sample_id of samples}
#' }
#' @param GTF The GTF file used to generate the peaks. This is used to determine the genomic coordinates of the gene.
#' @param ANNOTATION
#' @param OUTPUTDIR Output directory
#' @param OUTPUT.TAG A character string indicating a tag to track the generated files
#' @param PLOT Binary T or F indicating whether to save plot or just return plot
#'
#' @export plot.gene.peaks
#'
#' @examples
plot.gene.peaks = function(
  GENE,
  PEAKS,
  GTF = NULL,
  ANNOTATION = NULL,
  OUTPUTDIR = ".",
  OUTPUT.TAG = "",
  PLOT = T
){

  # Making a list of parameters to pass back and forth
  PARAMETERS = list()
  PARAMETERS$GENE = GENE
  PARAMETERS$GTF = GTF
  PARAMETERS$OUTPUTDIR = OUTPUTDIR
  PARAMETERS$OUTPUT.TAG = OUTPUT.TAG

  if(!PARAMETERS$GENE %in% PEAKS$name){
    warning("No Peaks are Found for This Gene in PEAKS!", call. = TRUE, domain = NULL)
    return(.generate.null.result(PARAMETERS))
  }

  # Import GTF as a GRanges Object
  annot.format = F
  if(!is.null(ANNOTATION)){
    PEAKSGR = .retrieve.peaks.as.granges(PEAKS = PEAKS, GENE = PARAMETERS$GENE, DF = F)
    annot.format = .check.annotation(ANNOTATION, PEAKSGR, GENE = PARAMETERS$GENE)
  }
  if(!annot.format){
    ANNOTATION = read.gtf(PARAMETERS)
  }

  # Plotting Peaks
  plotting.peaks = .retrieve.peaks.as.granges(PEAKS = PEAKS, GENE = PARAMETERS$GENE, DF = T)

  # Gene Bed
  gene.bed = ANNOTATION[ANNOTATION$gene == PARAMETERS$GENE, c("chr", "start", "stop", "gene", "strand")]
  gene.chr = unique(gene.bed$chr)

  # Code for ggplot
  options(warn = -1)
  p1 = ggplot2::ggplot(plotting.peaks, ggplot2::aes(y = sample, x = start, xend = end)) +
    ggalt::geom_dumbbell() +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 0),
                   axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::ggtitle(PARAMETERS$GENE) +
    ggplot2::xlab(gene.chr) + ggplot2::ylab("Sample")

  p1 = p1 + ggplot2::annotate(
    "rect",
    xmin=gene.bed$start,
    xmax=gene.bed$stop,
    ymin=-2,
    ymax=0,
    alpha=0.2,
    color="black",
    fill=rainbow(nrow(gene.bed))
    )
  options(warn = 0)

  if(PLOT){
    filename = paste0(PARAMETERS$OUTPUTDIR, "/", PARAMETERS$GENE, ".", PARAMETERS$OUTPUT.TAG, ".Peaks.pdf")
    pdf(filename)
    print(p1)
    dev.off()
  }

  return(p1)
}
