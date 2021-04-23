
#' simulate.gaussian.peaks
#'
#' @param MU An integer vector of means used to simulate a mixture of Gaussians. Values represent the mean bp, counted from the starting position of the transcript coordinate.
#' @param SD An integer vector of standard deviations used to simulate a mixture of Gaussians. Units are in bp.
#' @param EXTEND.WIDTH An integer vector of widths of peaks in bp.
#' @param NSAMPLES An integer vector
#' @param GENE The name of a gene in the GTF file on which to simulate peaks
#' @param GTF The path to a GTF file used to simulate RNA peaks
#' @param SEED A seed for reproducibility
#'
#' @return A data.frame in BED12 format, with an extra 'sample' column which can be used as input to ConsistentPeaks, PEAKS input
#'
#' @export simulate.gaussian.peaks
#'
#' @examples
#' gtf.path = system.file("extdata", package = "ConsensusPeaks")
#' gtf = paste0(gtf.path, "/test.gtf")
#' simulate.gaussian.peaks(
#' MU = c(100, 150),
#' SD = c(10, 20),
#' EXTEND.WIDTH = c(50, 25),
#' NSAMPLES = c(10, 30),
#' GENE = "ENSGXX",
#' GTF = gtf,
#' SEED = 123
#' )
#'
simulate.gaussian.peaks = function(
  MU = 1,
  SD = 1,
  EXTEND.WIDTH = 50,
  NSAMPLES = 10,
  GENE,
  GTF = NULL,
  ANNOTATION = NULL,
  SEED = 123
  ){

  # Error checking
  NGAUSSIANS = length(MU)
  params.lengths = c(length(MU), length(SD), length(EXTEND.WIDTH), length(NSAMPLES))
  if(sum(params.lengths) != 4*NGAUSSIANS){stop("MU, SD, EXTEND.WIDTH and NSAMPLES need to be the same length")}
  if(!all(MU %% 1 == 0)| !all(MU >= 0)){stop("MU must be an integer vector where all elements are greater than 0")}
  if(!all(SD %% 1 == 0) | !all(SD > 0)){stop("SD must be an integer vector where all elements are greater than 0")}
  if(!all(EXTEND.WIDTH %% 1 == 0) | !all(EXTEND.WIDTH > 0)){stop("EXTEND.WIDTH must be an integer vector where all elements are greater than 0")}
  if(!is.character(GENE) | GENE == "" | is.na(GENE)){stop("GENE must be character")}
  if(!(SEED %% 1 == 0)){stop("SEED must be integer")}

  if(is.null(ANNOTATION)) {
    # Creating annotation and geneinfo
    ANNOTATION = read.gtf(GTF)
  }
  if(!GENE %in% ANNOTATION$gene){stop("GENE must be in the GTF file")}
  GENEINFO = .get.gene.anno(GENE, ANNOTATION)

  # Initializing
  peaks = GenomicRanges::GRanges()
  set.seed(SEED)

  # Generating Peaks
  for(i in 1:NGAUSSIANS){

    # Creating RNA peaks
    pts = round(rnorm(NSAMPLES[i], mean = MU[i], sd = SD[i]))
    sample_ids = paste0("Sample.", i, ".", 1:NSAMPLES[i])
    rna.peaks = GenomicRanges::GRanges(seqnames = GENEINFO$chr, IRanges::IRanges(start = pts, end = pts), strand = GENEINFO$strand)
    GenomicRanges::mcols(rna.peaks)$sample = sample_ids
    GenomicRanges::mcols(rna.peaks)$name = rep(GENEINFO$gene, NSAMPLES[i])
    rna.peaks = GenomicRanges::resize(rna.peaks, width = EXTEND.WIDTH[i], fix = "center")
    GenomicRanges::start(rna.peaks) = ifelse(GenomicRanges::start(rna.peaks) < 1, 1, GenomicRanges::start(rna.peaks))
    GenomicRanges::end(rna.peaks) = ifelse(GenomicRanges::end(rna.peaks) > GENEINFO$exome_length, GENEINFO$exome_length, GenomicRanges::end(rna.peaks))
    peaks = c(peaks, rna.peaks)

  }

  # Creating DNA peaks
  dna.peaks = .rna.peaks.to.genome(peaks, GENEINFO)

  # Creating BED12 peaks
  dna.peaks.bed12  = .bed6tobed12(MERGED.PEAKS = dna.peaks, ID.COLS = "sample")
  colnames(dna.peaks.bed12)[which(colnames(dna.peaks.bed12) == "peak")] = "sample"

  return(dna.peaks.bed12)
}
