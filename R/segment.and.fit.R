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
  GENEGR = GenomicRanges::GRanges(seqnames = GENEINFO$chr, IRanges::IRanges(start = 1, end = GENEINFO$exome_length), strand = GENEINFO$strand)
  PEAK.COVERAGE = GenomicRanges::coverage(GENEPEAKSGR)
  BINS = unlist(GenomicRanges::tile(GENEGR, width = 1))
  BIN.COUNTS = data.frame(GenomicRanges::binnedAverage(BINS, PEAK.COVERAGE, "Coverage"), stringsAsFactors = F)

  # Segmenting & Determining which segments are peaks
  # Test 1: Fit smoothing spline before finding the peaks
  smooth.coverage = smooth.spline( BIN.COUNTS$start, BIN.COUNTS$Coverage, spar = 0.3)
  p = find.peaks(x = -smooth.coverage$y, m = 150, diff.threshold = 10^-7)
  # Test 2: As-is
  # p = find.peaks(-BIN.COUNTS$Coverage, m = 150)
  # Test 3: Moving Average
  # p.moving = moving.average(BIN.COUNTS$Coverage, n = 10)
  # p.moving[is.na(p.moving)] = BIN.COUNTS$Coverage[is.na(p.moving)]
  # p = find.peaks(-p.moving, m = 150)
  # Test 4: Remove Local Abnormalities
  # p.na.abnorm = remove.local.abnormalities(BIN.COUNTS, max_background_fold_increase = 2, background_window = 100)
  # p = find.peaks(-p.na.abnorm, m = 150)

  # Formatting
  SEG.GR = generate.peaks.from.split.points(p, GENEPEAKSGR, GENEINFO, m = 100)

  # Tiling Peaks
  peak.counts = unlist(GenomicRanges::tile(GENEPEAKSGR, width = 1))
  peak.counts = GenomicRanges::start(peak.counts)

  # Segments
  if(GENE %in% PARAMETERS$PLOT.MERGED.PEAKS) {
    filename = file.path(PARAMETERS$OUTPUTDIR, paste0(GENE, "segments.pdf"))
    pdf(filename, width = 5, height = 5)
    plot(BIN.COUNTS$start, BIN.COUNTS$Coverage, type = "s")
    lines(BIN.COUNTS$start, smooth.coverage$y, type = "s", col = "pink")
    # lines(BIN.COUNTS$start, p.moving, type = "s", col = "green")
    # lines(BIN.COUNTS$start, p.na.abnorm, type = "s", col = "chartreuse4")
    points(BIN.COUNTS$start[p], BIN.COUNTS$Coverage[p], col = 'red')
    dev.off()
  }

  # Fitting different models
  results = data.frame()
  models = list()
  for(i in 1:length(SEG.GR)){
    seg.start = GenomicRanges::start(SEG.GR)[i]
    seg.end = GenomicRanges::end(SEG.GR)[i]
    x = peak.counts[peak.counts >= seg.start & peak.counts <= seg.end]
    # Start the counts at 0 to fit the distributions better
    x.adj <- (x - seg.start) + 1e-10
    sd.scale <- sd(x.adj)
    x.scale <- x.adj / sd.scale
    mod = list()

    # Uniform Distribution
    tryCatch(
      expr = {
        mod$unif <- fitdistrplus::fitdist(
          data = x.scale,
          distr = "unif")
        mod$unif$sd_scale <- sd.scale
      },
      error = function(e) {
        warning(sprintf("Error in fitdist unif for segment [%d, %d]", seg.start, seg.end))
        print(e)
      }
    )

    # Attempt to fit a truncated normal distribution
    # Wrap in a tryCatch incase MLE fails to converge
    tryCatch(
      expr = {
        mod$tnorm <- fitdistrplus::fitdist(
          data = x.scale,
          distr = "tnorm",
          fix.arg = list(a = 0, b = max(x.scale) + 1e-10),
          start = list(mean = mean(x.scale), sd = sd(x.scale)),
          optim.method="L-BFGS-B")
        mod$tnorm$sd_scale <- sd.scale
      },
      error = function(e) {
        warning(sprintf("Error in fitdist tnorm for segment [%d, %d]", seg.start, seg.end))
        print(e)
      }
    )

    # Truncated gamma distribution
    tryCatch(
      expr = {
        mod$tgamma <- fitdistrplus::fitdist(
          data = x.scale,
          distr = "tgamma",
          fix.arg = list(a = 0, b = max(x.scale)),
          start = list(shape = 2, rate = 1),
          # Set the lower bound for shape and rate params
          lower = c(1, 0.5))
        mod$tgamma$sd_scale <- sd.scale
      },
      error = function(e) {
        warning(sprintf("Error in fitdist tgamma for segment [%d, %d]", seg.start, seg.end))
        print(e)
      }
    )

    # Flipped Truncated gamma distribution
    tryCatch(
      expr = {
        mod$tgamma_flip <- fitdistrplus::fitdist(
          data = x.scale,
          distr = "tgamma_flip",
          fix.arg = list(b = max(x.scale) + 1e-10),
          start = list(shape = 2, rate = 1),
          # Set the lower bound for shape and rate params
          lower = c(1, 0.5))
        mod$tgamma_flip$sd_scale <- sd.scale
      },
      error = function(e) {
        warning(sprintf("Error in fitdist tgamma_flip for segment [%d, %d]", seg.start, seg.end))
        print(e)
      }
    )

    fits = lapply(mod, function(mod){
      params <- c(mod$estimate, mod$fix.arg)

      data.frame(
        "dist" = summary(mod)$distname,
        "loglikelihood" = summary(mod)$loglik,
        "aic" = summary(mod)$aic,
        "bic" = summary(mod)$bic,
        # Text representation of the parameters
        params = dput.str(params)
      )
    })

    if(fit.norm_mixture) {
      # Fit a Normal Mixture model
      maxiter <- 500
      # Fit mixture model, silencing output
      invisible(capture.output({
        mixfit <- mixtools::normalmixEM(x.scale, verb = FALSE, maxit = maxiter, epsilon = 1e-04, k = 2)
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
          "bic" = -2 * mixfit$loglik + mixfit.params * log(length(x)),
          "params" = dput.str(mixfit[c("mu", "sigma", "lambda")])
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
  comput.fti <- function(i) { # models, SEG.GR
    # Initializing
    seg.start <- GenomicRanges::start(SEG.GR)[i]
    seg.end <- GenomicRanges::end(SEG.GR)[i]
    x = seg.start:seg.end
    fti = models[[i]]
    if(class(fti) == "mixEM") {
      scalefactor = length(fti$x)
      distname <- "norm_mixture"
      dens <- dnorm_mixture(x, fti)
    } else {
      # Center the region to match fitting process
      x.adj <- (x - seg.start) + 1e-10
      scalefactor <- length(fti$data)
      # Adjust the scale of the segment
      if(!is.null(fti$sd_scale)) {
        x.adj <- x.adj / fti$sd_scale
        # The bin size is now 1/sd_scale
        # Multiple scale factor by this new bin size
        scalefactor <- scalefactor / fti$sd_scale
      }

      # Generating data
      params <- c(as.list(fti$estimate), as.list(fti$fix.arg))
      distname <- fti$distname
      ddistname <- paste0("d", distname)
      call.params <- c(list(x = x.adj), as.list(params))
      dens <- do.call(ddistname, call.params)
    }
    data.frame(
      "x" = x,
      "dens" = dens * scalefactor,
      "col" = distname,
      row.names = NULL)
  }

  # Making a Nice Figure
  if(GENE %in% PARAMETERS$PLOT.MERGED.PEAKS) {

    distr.plotting.data = lapply(1:length(SEG.GR), comput.fti)

    filename = file.path(PARAMETERS$OUTPUTDIR, paste0(GENE, ".SegmentAndFit.pdf"))
    pdf(filename, width = 10, height = 10)

    p1 = ggplot2::ggplot(BIN.COUNTS, ggplot2::aes(x = start, y = Coverage)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(data = BIN.COUNTS[p,], ggplot2::aes(x = start, y = Coverage), col = 'red') +
      ggplot2::theme_bw() +
      ggplot2::ggtitle(GENEINFO$gene) +
      ggplot2::ylab("Coverage (at BP resolution)") +
      ggplot2::xlab("Transcript Coordinate") +
      ggplot2::annotate("rect", xmin=GenomicRanges::start(SEG.GR), xmax=GenomicRanges::end(SEG.GR), ymin=-1 , ymax=-0.1, alpha=0.5, color="black", fill=1:length(SEG.GR))

    for(i in 1:length(distr.plotting.data)){
      p1 = p1 + ggplot2::geom_line(data = distr.plotting.data[[i]], ggplot2::aes(x=x, y=dens, color = col))
    }
    p1 = p1 + ggplot2::guides(col=ggplot2::guide_legend(title="Distribution"))
    print(p1)

    dev.off()

  }

  # Generating Peaks
  merged.peaks = SEG.GR
  S4Vectors::mcols(merged.peaks)$i = results$i
  S4Vectors::mcols(merged.peaks)$dist = results$dist
  S4Vectors::mcols(merged.peaks)$params = results$params
  S4Vectors::mcols(merged.peaks)$name = GENEINFO$gene
  merged.peaks = .rna.peaks.to.genome(merged.peaks, GENEINFO)
  GenomicRanges::start(merged.peaks) = GenomicRanges::start(merged.peaks)-1

  # Generating BED12 File
  PEAKS.FINAL = .bed6tobed12(MERGED.PEAKS = merged.peaks, ID.COLS = c("name", "i", "dist"))
  # P-Value Table
  SAMPLE.PVAL = .merge.p(PEAKSGR, MERGED.PEAKS = merged.peaks, ANNOTATION, PARAMETERS, ID.COLS = c("name", "i", "dist"))
  # Output Table
  OUTPUT.TABLE = merge(PEAKS.FINAL, SAMPLE.PVAL, by = "peak", all = T)

  return(OUTPUT.TABLE)

}
