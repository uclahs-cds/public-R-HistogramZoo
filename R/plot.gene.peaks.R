
#' A plot of the peaks found for a particular gene
#'
#' @param gene A 'character' gene id corresponding to the gene_id's found in the gtf files and the 'name' column in the peak files
#' @param peaks A data frame containing the following columns, and potentially extras, usually found in a BED12 file, base 0 system
#' \describe{
#'   \item{chr}{chromosomes, character, same format as those identified in gtf file}
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
#' @param gtf The gtf file used to generate the peaks. This is used to determine the genomic coordinates of the gene.
#' @param output.dir Output directory
#' @param output.tag A character string indicating a tag to track the generated files
#' @param save.plot Binary T or F indicating whether to save plot or just return plot
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
#' gene = "ENSGXX",
#' peaks = peaks,
#' gtf = gtf,
#' output.dir = ".",
#' output.tag = "",
#' save.plot = F
#' )
plot.gene.peaks = function(
  gene,
  peaks,
  gtf = NULL,
  output.dir = ".",
  output.tag = "",
  save.plot = F
){

  if(!gene %in% peaks$name){stop("No Peaks are Found for This Gene in peaks!", call. = TRUE, domain = NULL)}

  # annotation
  annotation = read.gtf(gtf)

  # Plotting Peaks
  plotting.peaks = .retrieve.peaks.as.granges(peaks = peaks, gene = gene, DF = T)

  # Gene Bed
  gene.bed = annotation[annotation$gene == gene, c("chr", "start", "stop", "gene", "strand")]
  gene.chr = unique(gene.bed$chr)

  # Code for ggplot
  options(warn = -1)
  p1 = ggplot2::ggplot(plotting.peaks, ggplot2::aes(y = sample, x = start, xend = end)) +
    ggalt::geom_dumbbell() +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 0),
                   axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::ggtitle(gene) +
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

  if(save.plot){
    filename = file.path(output.dir, paste0(gene, ".", output.tag, ".Peaks.pdf"))
    pdf(filename)
    print(p1)
    dev.off()
  }

  return(p1)
}
