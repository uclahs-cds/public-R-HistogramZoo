
#' A plot of the peaks found for a particular gene
#'
#' @param GENE A 'character' gene id corresponding to the gene_id's found in the GTF files and the 'name' column in the peak files
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
#' @param GTF The GTF file used to generate the peaks. This is used to determine the genomic coordinates of the gene.
#' @param OUTPUTDIR Output directory
#' @param OUTPUT.TAG A character string indicating a tag to track the generated files
#' @param PLOT Binary T or F indicating whether to save plot or just return plot
#'
#' @export plot.gene.peaks
#'
#' @return A ggplot object of the plot
#'
#' @examples
#' data.path = system.file("extdata", package = "ConsensusPeaks")
#' gtf = paste0(data.path, "/test.gtf")
#' peaks.file = paste0(data.path, "/peaks.bed")
#' peaks = read.delim(peaks.file, stringsAsFactors = F)
#' p = plot.gene.peaks(
#' GENE = "ENSGXX",
#' PEAKS = peaks,
#' GTF = gtf,
#' OUTPUTDIR = ".",
#' OUTPUT.TAG = "",
#' PLOT = F
#' )
plot.gene.peaks = function(
  GENE,
  PEAKS,
  GTF = NULL,
  OUTPUTDIR = ".",
  OUTPUT.TAG = "",
  PLOT = F
){

  if(!GENE %in% PEAKS$name){stop("No Peaks are Found for This Gene in PEAKS!", call. = TRUE, domain = NULL)}

  # ANNOTATION
  ANNOTATION = read.gtf(GTF)

  # Plotting Peaks
  plotting.peaks = .retrieve.peaks.as.granges(PEAKS = PEAKS, GENE = GENE, DF = T)

  # Gene Bed
  gene.bed = ANNOTATION[ANNOTATION$gene == GENE, c("chr", "start", "stop", "gene", "strand")]
  gene.chr = unique(gene.bed$chr)

  # Code for ggplot
  options(warn = -1)
  p1 = ggplot2::ggplot(plotting.peaks, ggplot2::aes(y = sample, x = start, xend = end)) +
    ggalt::geom_dumbbell() +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 0),
                   axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::ggtitle(GENE) +
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
    filename = paste0(OUTPUTDIR, "/", GENE, ".", OUTPUT.TAG, ".Peaks.pdf")
    pdf(filename)
    print(p1)
    dev.off()
  }

  return(p1)
}
