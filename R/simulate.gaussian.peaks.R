
#' simulate.gaussian.peaks
#'
#' @param mu An integer vector of means used to simulate a mixture of Gaussians. Values represent the mean bp, counted from the starting position of the transcript coordinate.
#' @param sd An integer vector of standard deviations used to simulate a mixture of Gaussians. Units are in bp.
#' @param extend.width An integer vector of widths of peaks in bp.
#' @param nsamples An integer vector
#' @param gene The name of a gene in the gtf file on which to simulate peaks
#' @param gtf The path to a gtf file used to simulate RNA peaks
#' @param seed A seed for reproducibility
#'
#' @return A data.frame in BED12 format, with an extra 'sample' column which can be used as input to ConsistentPeaks, PEAKS input
#'
#' @export simulate.gaussian.peaks
#'
#' @examples
#' gtf.path = system.file("extdata", package = "ConsensusPeaks")
#' gtf = paste0(gtf.path, "/test.gtf")
#' simulate.gaussian.peaks(
#' mu = c(100, 150),
#' sd = c(10, 20),
#' extend.width = c(50, 25),
#' nsamples = c(10, 30),
#' gene = "ENSGXX",
#' gtf = gtf,
#' seed = 123
#' )
#'
simulate.gaussian.peaks = function(
  mu = 1,
  sd = 1,
  extend.width = 50,
  nsamples = 10,
  gene,
  gtf = NULL,
  ANNOTATION = NULL,
  seed = 123
  ){

  # Error checking
  NGAUSSIANS = length(mu)
  params.lengths = c(length(mu), length(sd), length(extend.width), length(nsamples))
  if(sum(params.lengths) != 4*NGAUSSIANS){stop("mu, sd, extend.width and nsamples need to be the same length")}
  if(!all(mu %% 1 == 0)| !all(mu >= 0)){stop("mu must be an integer vector where all elements are greater than 0")}
  if(!all(sd %% 1 == 0) | !all(sd > 0)){stop("sd must be an integer vector where all elements are greater than 0")}
  if(!all(extend.width %% 1 == 0) | !all(extend.width > 0)){stop("extend.width must be an integer vector where all elements are greater than 0")}
  if(!is.character(gene) | gene == "" | is.na(gene)){stop("gene must be character")}
  if(!(seed %% 1 == 0)){stop("seed must be integer")}

  if(is.null(ANNOTATION)) {
    # Creating annotation and geneinfo
    ANNOTATION = read.gtf(gtf)
  }
  if(!gene %in% ANNOTATION$gene){stop("gene must be in the gtf file")}
  geneINFO = .get.gene.anno(gene, ANNOTATION)

  # Initializing
  peaks = GenomicRanges::GRanges()
  set.seed(seed)

  # Generating Peaks
  for(i in 1:NGAUSSIANS){

    # Creating RNA peaks
    pts = round(rnorm(nsamples[i], mean = mu[i], sd = sd[i]))
    sample_ids = paste0("Sample.", i, ".", 1:nsamples[i])
    rna.peaks = GenomicRanges::GRanges(seqnames = geneINFO$chr, IRanges::IRanges(start = pts, end = pts), strand = geneINFO$strand)
    GenomicRanges::mcols(rna.peaks)$sample = sample_ids
    GenomicRanges::mcols(rna.peaks)$name = rep(geneINFO$gene, nsamples[i])
    rna.peaks = GenomicRanges::resize(rna.peaks, width = extend.width[i], fix = "center")
    GenomicRanges::start(rna.peaks) = ifelse(GenomicRanges::start(rna.peaks) < 1, 1, GenomicRanges::start(rna.peaks))
    GenomicRanges::end(rna.peaks) = ifelse(GenomicRanges::end(rna.peaks) > geneINFO$exome_length, geneINFO$exome_length, GenomicRanges::end(rna.peaks))
    peaks = c(peaks, rna.peaks)

  }

  # Creating DNA peaks
  dna.peaks = .rna.peaks.to.genome(peaks, geneINFO)

  # Creating BED12 peaks
  dna.peaks.bed12  = .bed6tobed12(MERGED.PEAKS = dna.peaks, ID.COLS = "sample")
  colnames(dna.peaks.bed12)[which(colnames(dna.peaks.bed12) == "peak")] = "sample"

  return(dna.peaks.bed12)
}
