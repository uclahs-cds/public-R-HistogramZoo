.get.gene.anno <- function(GENE, ANNOTATION) {

  # extract batch annotation
  anno=ANNOTATION[ANNOTATION$gene == GENE,]
  anno_unique=unique(anno)

  # extract information
  strand=as.character(anno_unique[1, "strand"])
  chr=as.character(anno_unique[1, "chr"])
  left=min(anno_unique["start"])
  right=max(anno_unique["stop"])
  intervals=anno_unique[,c("start", "stop")]-left+1
  dna_length=right-left+1

  # prepare DNA2RNA
  DNA2RNA=rep(0,dna_length)
  no_intervals=length(intervals[,1])
  for (i in 1:no_intervals) {DNA2RNA[intervals[i,1]:intervals[i,2]]=1}
  exome_length=sum(DNA2RNA) # this is actually exome length
  DNA2RNA=cumsum(DNA2RNA)*DNA2RNA

  # prepare RNA2DNA
  RNA2DNA = left:right
  RNA2DNA = RNA2DNA[DNA2RNA > 0]

  # summarize result
  batch_anno=list(anno=anno_unique, gene=GENE,chr=chr,strand=strand,left=left,right=right,
                  DNA2RNA=DNA2RNA, RNA2DNA=RNA2DNA, dna_length=dna_length,
                  exome_length=exome_length)
}
