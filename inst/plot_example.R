#' Returns the group matches from a regular expression on a vector
#' @param x the vector we want to match on
#' @param pattern regular expression with groups
str.match <- function(x, pattern) {
  regmatches(x, regexec(pattern, x));
}

read_peak <- function(filename, sample) {
  if(missing(sample)) {
    sample <- str.match(filename, "CPCG[0-9]{4}")[[1]][1]
  }
  peak_df <- read.table(filename, header = F, sep = "\t", stringsAsFactors = F)
  peak_df$sample <- sample
  colnames(peak_df) = c("chr", "start", "end", "name", "score", "strand", "thickStart", "thickEnd",
                        "itemRGB", "blockCount", "blockSizes", "blockStarts", "sample")
  peak_df
}

gtf.file <- '~/data/PeakFitting/gencode.v34.chr_patch_hapl_scaff.annotation.gtf'

PARAMETERS <- list()
PARAMETERS$GTF <- gtf.file
annotation <- read.gtf(PARAMETERS)

met_peak_files <- list.files("~/data/PeakFitting/F_MeTPeak_Peaks_SS/", full.names = T, recursive = T)
exome_peak_files <- list.files("~/data/PeakFitting/G_exomePeak_Peaks_SS/", full.names = T, recursive = T)

gene <- "ENSG00000142515.15"

met_peaks <- lapply(met_peak_files, read_peak)
met_peaks <- do.call(rbind.data.frame, met_peaks)

exome_peaks <- lapply(exome_peak_files, read_peak)
exome_peaks <- do.call(rbind.data.frame, exome_peaks)

gene.coverage.hist(exome_peaks, gene, annotation) + ggtitle(paste("Exome Peaks", gene))
gene.coverage.hist(met_peaks, gene, annotation) + ggtitle(paste("MeT Peaks", gene))
