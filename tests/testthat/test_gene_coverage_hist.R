context("Testing Plot Gene Coverage")

gtf.path = system.file("extdata", package = "ConsensusPeaks")
gtf = paste0(gtf.path, "/test.gtf")
PARAMETERS = list("GTF" = gtf, "GENE" = "ENSGXX")
annotation = read.gtf(PARAMETERS)

ensgxx = simulate.gaussian.peaks(
  MU = c(100, 150),
  SD = c(10, 20),
  EXTEND.WIDTH = c(50, 25),
  NSAMPLES = c(10, 30),
  GENE = "ENSGXX",
  GTF = gtf,
  ANNOTATION = annotation,
  SEED = 123)

gene.coverage.hist(ensgxx, "ENSGXX", annotation)
