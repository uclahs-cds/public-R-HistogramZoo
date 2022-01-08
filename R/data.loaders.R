

# A function to produce a gene model out of a GTF file
gtf.to.genemodel = function(
  gtf.file,
  gene.or.transcript = c("gene", "transcript"),
  select.chrs = NULL
){
  gtf = valr::read_gtf(gtf.file)
  gtf = gtf[gtf$type == "exon" & gtf$chrom %in% select.chrs,]
  if(gene.or.transcript == "gene"){
    gtf[,"id"] = gtf[,"gene_id"]
  } else {
    gtf[,"id"] = gtf[,"transcript_id"]
  }
  gtf = gtf[,c("chrom", "start", "end", "strand", "id")]
  gtf.gr = GenomicRanges::makeGRangesFromDataFrame(
    gtf, 
    keep.extra.columns = T,
    seqinfo = select.chrs
  )
  gtf.gr = S4Vectors::split(gtf.gr, gtf.gr$id)
  gtf.gr = GenomicRanges::reduce(gtf.gr)
  gtf.gr
}

# A function that takes a coverage Rle object and a GRangesList
# and generates a set of histograms
# The following part can be written into a need-to-use basis
# It might be more efficient than the current set-up
# *** Ask Stefan *** He probably knows
coverage.to.hist = function(
  regions,
  cov,
  histogram.bin.size
){
  histogram.coverage = list()
  for(i in 1:length(regions)){
    x = regions[i]
    x.name = names(x)
    x = unlist(x)
    bins = GenomicRanges::tile(x = x, width = histogram.bin.size)
    bins = unlist(bins)
    cvg = GenomicRanges::binnedAverage(
      bins = bins, 
      numvar = cov, 
      varname = "cvg")
    cvg = cvg$cvg
    names(cvg) = NULL
    histogram.coverage[[x.name]] <- cvg
  }
  histogram.coverage
}

# Imports BED files and exports something usable by the segmentation function
bed.to.hist = function(
  filenames, 
  n_fields = 3, 
  gtf.file = NULL,
  gene.or.transcript = c("gene", "transcript"),
  histogram.bin.size = 1
){
  
  # Load BED files
  bedfiles = valr::read_bed( filenames, n_fields = n_fields )
  if(n_fields == 12){ segs = valr::bed12_to_exons( bedfiles ) }
  segs.gr = GenomicRanges::makeGRangesFromDataFrame( segs )
  
  # Calculating coverage using GenomicRanges
  segs.coverage = GenomicRanges::coverage( segs.gr )
  
  # GTF file for RNA-based coordinates
  if(!is.null(gtf.file)){
    regions.gr = gtf.to.genemodel(
      gtf.file = gtf.file, 
      gene.or.transcript = gene.or.transcript,
      select.chrs = names(segs.coverage))
   
    # Filtering regions.gr 
    select.regions = S4Vectors::subjectHits(GenomicRanges::findOverlaps(segs.gr, regions.gr))
    regions.gr = regions.gr[sort(unique(select.regions))]
  } else {
    # Genome coordinates
    # Rethink this a little bit***
    # regions.gr = GenomicRanges::reduce(segs.gr)
    # regions.gr = GenomicRanges::split(regions.gr, 1:length(regions.gr))
  }
  
  histogram.coverage = coverage.to.hist(
    cov = segs.coverage, 
    regions = regions.gr, 
    histogram.bin.size = histogram.bin.size)
  
  # Returning list of histograms, region coordinates
  list(
    histogram.coverage = histogram.coverage, 
    gene.model = regions.gr, 
    histogram.bin.size = histogram.bin.size)
}

# Imports BigWig files and exports something usable by the segmentation function
bigwig.to.hist = function(
  filename,
  set.strand = c("+", "-"),
  gtf.file = NULL,
  gene.or.transcript = c("gene", "transcript"),
  histogram.bin.size = 1
){
  
  # Load BigWig file
  bigwig = valr::read_bigwig(filename, set_strand = set.strand)
  bigwig.gr = GenomicRanges::makeGRangesFromDataFrame(bigwig, keep.extra.columns = T)
  segs.coverage = GenomicRanges::coverage(bigwig.gr, weight = "score")
  bigwig.gr = bigwig.gr[bigwig.gr$score != 0]
  
  # Loading regions
  if(!is.null(gtf.file)){
    regions.gr = gtf.to.genemodel(
      gtf.file = gtf.file, 
      gene.or.transcript = gene.or.transcript,
      select.chrs = names(segs.coverage))
    
    # Filtering regions.gr 
    select.regions = S4Vectors::subjectHits(GenomicRanges::findOverlaps(bigwig.gr, regions.gr))
    regions.gr = regions.gr[sort(unique(select.regions))]
  } else {
    # Genome coordinates
    # Rethink this a little bit***
    # regions.gr = GenomicRanges::reduce(cov.gr)
    # regions.gr = GenomicRanges::split(regions.gr, 1:length(regions.gr))
  }
  
  histogram.coverage = coverage.to.hist(
    cov = segs.coverage, 
    regions = regions.gr, 
    histogram.bin.size = histogram.bin.size)
  
  # Returning list of histograms, region coordinates
  list(
    histogram.coverage = histogram.coverage, 
    gene.model = regions.gr, 
    histogram.bin.size = histogram.bin.size)
}


# Testing -----------------------------------------------------------------

library(valr)
library(genomation)
library(GenomicRanges)

setwd("/cluster/home/helenzhu/Cluster_Helen/Snakemake_ConsensusPeaks/TestData")
filenames = c("s1.bed", "s2.bed", "s3.bed")
n_fields = 12
gtf.file = "gencode.v34.chr_patch_hapl_scaff.annotation.gtf"
gene.or.transcript = "gene"
filename = "S1.bw"
set.strand = "+"

res = bed.to.hist(
  filenames = filenames,
  n_fields = 12,
  gtf.file = gtf.file,
  gene.or.transcript = "gene",
  histogram.bin.size = 1
)

res2 = bigwig.to.hist(
  filename = "S1.bw",
  set.strand = "+",
  gtf.file = gtf.file,
  gene.or.transcript = "gene",
  histogram.bin.size = 1
)

