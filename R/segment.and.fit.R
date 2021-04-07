#' Segment and fit a gene
#'
#' @param GENE
#'
#' @param PARAMETERS
#' @param ANNOTATION
#' @param PEAKS
#'
#' @import extraDistr
#' @export
segment.and.fit = function(
  GENE,
  PARAMETERS,
  ANNOTATION,
  PEAKS
){

  # If the gene doesn't have peaks
  if(!GENE %in% PEAKS$name){
    warn.message = paste0("No Peaks are Found for ", GENE, " in PEAKS!")
    warning(warn.message, call. = TRUE, domain = NULL)
    return(.generate.null.result(PARAMETERS))
  }

  # Use the parameters if they are defined and default to FALSE if not defined
  # plot.merged.peak <- PARAMETERS$PLOT.MERGED.PEAKS  %||% FALSE
  plot.diagnostic <- PARAMETERS$DIAGNOSTIC %||% FALSE
  write.output <- PARAMETERS$WRITE.OUTPUT %||% FALSE
  # Optional normal mixture model via mixtools
  fit.norm_mixture <- PARAMETERS$FIT.MIXTURE %||% FALSE

  # PEAKSGR
  PEAKSGR = ConsensusPeaks:::.retrieve.peaks.as.granges(PEAKS = PEAKS, GENE = GENE, DF = F)

  # Get Gene Information
  GENEINFO = ConsensusPeaks:::.get.gene.anno(GENE, ANNOTATION)

  # Converting to RNA
  GENEPEAKSGR = GenomicRanges::shift(PEAKSGR, -1*GENEINFO$left+1)
  GenomicRanges::start(GENEPEAKSGR) = GENEINFO$DNA2RNA[GenomicRanges::start(GENEPEAKSGR)+1]
  GenomicRanges::end(GENEPEAKSGR) = GENEINFO$DNA2RNA[GenomicRanges::end(GENEPEAKSGR)]

  # Reduce Overlapping Peaks in the Same Sample
  GENEPEAKSGR = S4Vectors::split(GENEPEAKSGR, GENEPEAKSGR$sample)
  GENEPEAKSGR = unlist(GenomicRanges::reduce(GENEPEAKSGR))

  # Peak Coverage
  REDUCED.GENE.PEAKS.GR = GenomicRanges::reduce(GENEPEAKSGR)
  PEAK.COVERAGE = GenomicRanges::coverage(GENEPEAKSGR)
  BINS = unlist(GenomicRanges::tile(REDUCED.GENE.PEAKS.GR, width = 1))
  BIN.COUNTS = data.frame(GenomicRanges::binnedAverage(BINS, PEAK.COVERAGE, "Coverage"), stringsAsFactors = F)
  add.zero.counts = setdiff(1:GENEINFO$exome_length, BIN.COUNTS$start)
  add.zero.counts = data.frame("start" = add.zero.counts, "Coverage" = 0, stringsAsFactors = F)
  BIN.COUNTS = rbind(BIN.COUNTS[,c("start", "Coverage")], add.zero.counts)
  BIN.COUNTS = BIN.COUNTS[order(BIN.COUNTS$start),]

  # Segmenting & Determining which segments are peaks
  # Test 1: Fit smoothing spline before finding the peaks
  # smooth.coverage = smooth.spline( BIN.COUNTS$start, BIN.COUNTS$Coverage, spar = 0.5)
  # p = find.peaks(-smooth.coverage$y, m = 150)
  # Test 2: As-is
  p = find.peaks(-BIN.COUNTS$Coverage, m = 150)
  # Test 3: Moving Average
  # p.moving = moving.average(BIN.COUNTS$Coverage, n = 10)
  # p = find.peaks(-p.moving, m = 150)

  # Formatting
  p = c(1, GENEINFO$exome_length, p)
  p = unique(p)
  p = sort(p)
  # Remove segments that are less than 100 apart
  p = p[diff(p) > 100]
  seg.df = data.frame("start" = p[1:(length(p)-1)]+1, "end" = p[2:length(p)], "mean" = 0)
  for(i in 1:nrow(seg.df)){
    tmp = BIN.COUNTS$Coverage[BIN.COUNTS$start >= seg.df$start[i] & BIN.COUNTS$start <= seg.df$end[i]]
    seg.df$mean[i] = mean(tmp)
  }
  filter.cond <- seg.df$mean > 0
  seg.df = seg.df[filter.cond, ]

  # Tiling Peaks
  peak.counts = unlist(GenomicRanges::tile(GENEPEAKSGR, width = 1))
  peak.counts = GenomicRanges::start(peak.counts)

  # Segments
  if(GENE %in% PARAMETERS$PLOT.MERGED.PEAKS) {
    filename = file.path(PARAMETERS$OUTPUTDIR, paste0(GENE, "segments.pdf"))
    pdf(filename, width = 5, height = 5)
    plot(BIN.COUNTS$start, BIN.COUNTS$Coverage, type = "s")
    # lines(BIN.COUNTS$start, smooth.coverage$y, type = "s", col = "pink")
    # lines(BIN.COUNTS$start, p.moving, type = "s", col = "green")
    points(BIN.COUNTS$start[p], BIN.COUNTS$Coverage[p], col = 'red')
    dev.off()
  }

  # Fitting different models
  results = data.frame()
  models = list()
  for(i in 1:nrow(seg.df)){
    x = peak.counts[peak.counts >= seg.df$start[i] & peak.counts <= seg.df$end[i]]
    simple.model.names <- c("norm", "gamma", "unif")

    mod = lapply(simple.model.names, function(distr) fitdistrplus::fitdist(x, distr))
    names(mod) <- simple.model.names

    # Attempt to fit a truncated normal distribution
    # Wrap in a tryCatch incase MLE fails to converge
    tryCatch(
      expr = {
        mod$tnorm <- fitdistrplus::fitdist(
          data = x,
          distr = "tnorm",
          fix.arg = list(a = seg.df$start[i] - 1, b = seg.df$end[i] + 1),
          start = list(mean = mean(x), sd = sd(x)),
          optim.method="L-BFGS-B")
      },
      error = function(e) {
        warning(sprintf("Error in fitdist tnorm for segment [%d, %d]", seg.df$start[i], seg.df$end[i]))
        print(e)
      }
    )

    fits = lapply(mod, function(mod){
      data.frame(
        "dist" = summary(mod)$distname,
        "loglikelihood" = summary(mod)$loglik,
        "aic" = summary(mod)$aic,
        "bic" = summary(mod)$bic
      )
    })

    if(fit.norm_mixture) {
      # Fit a Normal Mixture model
      maxiter <- 500
      # Fit mixture model, silencing output
      invisible(capture.output({
        mixfit <- mixtools::normalmixEM(x, verb = FALSE, maxit = maxiter, epsilon = 1e-04)
      }))
      if((length(mixfit$all.loglik) - 1) >= maxiter) {
        # EM did not converge. Don't use results.
        mixfit <- NULL
      } else {
        mixfit.params <- length(mixfit$mu) * 3 # for lambda, mu, and sigma params
        mixfit.results <- data.frame(
          "dist" = "norm_mixture",
          "loglikelihood" = mixfit$loglik,
          "aic" = -2 * mixfit$loglik + 2 * mixfit.params,
          "bic" = -2 * mixfit$loglik + mixfit.params * log(length(x))
        )
        fits$norm_mixture <- mixfit.results
      }
    }

    fits = do.call(rbind, fits)
    fits$i = i

    if(plot.diagnostic) {
      filename = file.path(PARAMETERS$OUTPUTDIR, paste0(GENE, ".fit.segments.", i, ".pdf"))
      pdf(filename, width = 10, height = 10)

      par(mfrow = c(2,2))
      plot.legend = names(mod)

      fitdistrplus::denscomp(mod, legendtext = plot.legend)
      fitdistrplus::qqcomp(mod, legendtext = plot.legend)
      fitdistrplus::cdfcomp(mod, legendtext = plot.legend)
      fitdistrplus::ppcomp(mod, legendtext = plot.legend)
      dev.off()
    }

    # Add the normal mixture after fitdistrplus
    if(fit.norm_mixture && !is.null(mixfit)) mod$norm_mixture <- mixfit
    # Adding the results to the table
    res.final = fits[fits$aic == min(fits$aic),]
    results = rbind(results, res.final)
    mod.final = mod[[res.final$dist]]
    models[[i]] = mod.final
  }

  # Acquiring a Density Distribution
  comput.fti <- function(i) { # models, seg.df
    # Initializing
    x = seg.df$start[i]:seg.df$end[i]
    fti = models[[i]]
    if(class(fti) == "mixEM") {
      scalefactor = length(fti$x)
      distname <- "norm_mixture"
      dens <- dnorm_mixture(x, fti)
    } else {
      scalefactor = length(fti$data)
      # Scale factor
      # Generating data
      para <- c(as.list(fti$estimate), as.list(fti$fix.arg))
      distname <- fti$distname
      ddistname <- paste("d", distname, sep = "")
      dens <- do.call(ddistname, c(list(x), as.list(para)))
    }
    data.frame("x" = x, "dens" = dens * scalefactor, "col" = distname, row.names = NULL)
  }

  # Making a Nice Figure
  if(GENE %in% PARAMETERS$PLOT.MERGED.PEAKS) {

    distr.plotting.data = lapply(1:nrow(seg.df), comput.fti)

    filename = file.path(PARAMETERS$OUTPUTDIR, paste0(GENE, ".SegmentAndFit.pdf"))
    pdf(filename, width = 10, height = 10)

    p1 = ggplot2::ggplot(BIN.COUNTS, ggplot2::aes(x = start, y = Coverage)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(data = BIN.COUNTS[p,], ggplot2::aes(x = start, y = Coverage), col = 'red') +
      ggplot2::theme_bw() +
      ggplot2::ggtitle(GENEINFO$gene) +
      ggplot2::ylab("Coverage (at BP resolution)") +
      ggplot2::xlab("Transcript Coordinate") +
      ggplot2::annotate("rect", xmin=seg.df$start, xmax=seg.df$end, ymin=-1 , ymax=-0.1, alpha=0.5, color="black", fill=1:nrow(seg.df))

    for(i in 1:length(distr.plotting.data)){
      p1 = p1 + ggplot2::geom_line(data = distr.plotting.data[[i]], ggplot2::aes(x=x, y=dens, color = col))
    }
    p1 = p1 + ggplot2::guides(col=ggplot2::guide_legend(title="Distribution"))
    print(p1)

    dev.off()

  }

  # Generating Peaks
  merged.peaks = GenomicRanges::GRanges(seqnames = GENEINFO$chr, IRanges::IRanges(start = seg.df$start, end = seg.df$end), strand = GENEINFO$strand, i = results$i, dist = results$dist, name = GENEINFO$gene)
  merged.peaks = ConsensusPeaks:::.rna.peaks.to.genome(merged.peaks, GENEINFO)
  GenomicRanges::start(merged.peaks) = GenomicRanges::start(merged.peaks)-1

  # Generating BED12 File
  PEAKS.FINAL = ConsensusPeaks:::.bed6tobed12(MERGED.PEAKS = merged.peaks, ID.COLS = c("name", "i", "dist"))
  # P-Value Table
  SAMPLE.PVAL = ConsensusPeaks:::.merge.p(PEAKSGR, MERGED.PEAKS = merged.peaks, ANNOTATION, PARAMETERS, ID.COLS = c("name", "i", "dist"))
  # Output Table
  OUTPUT.TABLE = merge(PEAKS.FINAL, SAMPLE.PVAL, by = "peak", all = T)

  return(OUTPUT.TABLE)

}
