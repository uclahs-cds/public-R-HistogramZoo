
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
#' @export
#'
#' @examples
#' gtf.path = system.file("extdata", package = "DPDE4PM")
#' gtf = paste0(gtf.path, "/test.gtf")
#' simulate.gaussian.peaks(
#' MU = c(100, 150),
#' SD = c(10, 20),
#' EXTEND.WIDTH = c(50, 25),
#' NSAMPLES = c(10, 30),
#' GENE = "",
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
  SEED = 123
  ){
  
  # Error checking
  NGAUSSIANS = length(MU)
  params.lengths = c(length(MU), length(SD), length(EXTEND.WIDTH), length(NSAMPLES))
  if(sum(params.lengths) != 4*NGAUSSIANS){stop("MU, SD, EXTEND.WIDTH and NSAMPLES need to be the same length")}
  if(!is.integer(MU) | !all(MU >= 0)){stop("MU must be an integer vector where all elements are greater than 0")}
  if(!is.integer(SD) | !all(SD > 0)){stop("SD must be an integer vector where all elements are greater than 0")}
  if(!is.integer(EXTEND.WIDTH) | !all(EXTEND.WIDTH > 0)){stop("EXTEND.WIDTH must be an integer vector where all elements are greater than 0")}
  if(!is.character(GENE)){stop("GENE must be character")}
  if(!is.integer(SEED)){stop("SEED must be integer")}
  
  # PARAMETERS 
  PARAMETERS = list()
  PARAMETERS$GENE = GENE
  PARAMETERS$GTF = GTF
  
  # Creating annotation and geneinfo
  ANNOTATION = read.gtf(PARAMETERS)
  GENEINFO = .get.gene.anno(PARAMETERS, ANNOTATION)
  
  # Initializing
  peak = GenomicRanges::GRanges()
  set.seed(SEED)
  
  # Generating Peaks
  for(i in 1:NGAUSSIANS){
    
    # Creating RNA peaks
    pts = round(rnorm(NSAMPLES[i], mean = MU[i], sd = SD[i]))
    sample_ids = generate.random.ids(NSAMPLES[i])
    rna.peaks = GenomicRanges::GRanges(seqnames = GENEINFO$chr, IRanges::IRanges(start = pts, end = pts), mcols = sample_ids, strand = GENEINFO$strand)
    rna.peaks = GenomicRanges::resize(rna.peaks, width = EXTEND.WIDTH[i], fix = "center")
    GenomicRanges::start(rna.peaks) = ifelse(GenomicRanges::start(rna.peaks) < 1, 1, GenomicRanges::start(rna.peaks))
    GenomicRanges::end(rna.peaks) = ifelse(GenomicRanges::end(rna.peaks) > GENEINFO$exome_length, GENEINFO$exome_length, GenomicRanges::end(rna.peaks))
    peaks = c(peaks, rna.peaks)
    
  }
  
  # Creating DNA peaks
  dna.peaks = .rna.peaks.to.genome(peaks, GENEINFO)
  
  # Creating BED12 peaks
  dna.peaks.bed12  = .bed6tobed12(MERGED.PEAKS = dna.peaks, ID.COLS = "mcols")
  
  return(dna.peaks.bed12)
}

