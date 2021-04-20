#' read.gtf is a function that imports a gtf file and creates an annotation data frame, adapted from exomePeak
#'
#' @param PARAMETERS A PARAMETERS list with the parameters indicated in the DPDE4PM function
#'
#' @return
#' \describe{
#'  \item{chr}{chromosome}
#'  \item{feature}{genomic feature}
#'  \item{start}{start coordinate, base 1}
#'  \item{end}{stop coordinate, base 1}
#'  \item{strand}{strand}
#'  \item{gene}{gene id in the GTF file}
#'  \item{transcript}{transcript id in the GTF fiile}
#' }
#'
#' @export read.gtf
read.gtf <- function(PARAMETERS){

  # Creating a TXDB
  txdb=suppressWarnings(GenomicFeatures::makeTxDbFromGFF(PARAMETERS$GTF,format="gtf"))

  # Filtering the TXDB
  colkey <- AnnotationDbi::columns(txdb)
  select_col <- match(c("EXONCHROM","TXID","EXONSTART","EXONEND","EXONSTRAND","GENEID","TXNAME"),colkey)
  ID = AnnotationDbi::keys(txdb, "TXID")
  temp = AnnotationDbi::select(txdb, ID , c(AnnotationDbi::columns(txdb))[select_col], "TXID")
  select_col2 <- match(c("EXONCHROM","TXID","EXONSTART","EXONEND","EXONSTRAND","GENEID","TXNAME"),names(temp))
  temp <- temp[,select_col2]
  colnames(temp)=c("chr","feature","start","stop","strand","gene","transcript")
  temp$"feature" <- "exon"
  gtf <- temp

  # return data
  return(gtf)
}
