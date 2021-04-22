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
  PEAKSGR = .retrieve.peaks.as.granges(PEAKS = PEAKS, GENE = GENE, DF = F)

  # Get Gene Information
  GENEINFO = .get.gene.anno(GENE, ANNOTATION)

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
  smooth.coverage = smooth.spline( BIN.COUNTS$start, BIN.COUNTS$Coverage, spar = 0.3)
  p = find.peaks(x = -smooth.coverage$y, m = 150, diff.threshold = 10^-7)

  # Formatting
  SEG.GR = generate.peaks.from.split.points(p, GENEPEAKSGR, GENEINFO, m = 100)

  # Tiling Peaks
  peak.counts = unlist(GenomicRanges::tile(GENEPEAKSGR, width = 1))
  peak.counts = GenomicRanges::start(peak.counts)

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

    x.range = seg.start:seg.end
    x.range.adj <- (x - seg.start) + 1e-10
    # The bin size is now 1/sd_scale
    x.range.scale <- x.range.adj / sd.scale
    # Multiple scale factor by this new bin size
    scalefactor <- length(x.scale) / sd.scale

    fits = lapply(mod, function(m){
      params <- c(m$estimate, m$fix.arg)

      bin.data <- table(x.scale)

      params <- c(as.list(m$estimate), as.list(m$fix.arg))
      distname <- m$distname
      ddistname <- paste0("d", distname)
      call.params <- c(list(x = as.numeric(names(bin.data))), as.list(params))
      dens <- do.call(ddistname, call.params)
      dens.scale <- dens * scalefactor

      fit.residuals <- (dens.scale - as.integer(bin.data))^2

      data.frame(
        "dist" = summary(m)$distname,
        "loglikelihood" = summary(m)$loglik,
        "aic" = summary(m)$aic,
        "bic" = summary(m)$bic,
        "mse" = mean(fit.residuals),
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

        bin.data <- table(x.scale)
        dens <- dnorm_mixture(as.numeric(names(bin.data)), mixfit)
        dens.scale <- dens * scalefactor

        fit.residuals <- (dens.scale - as.integer(bin.data))^2

        mixfit.results <- data.frame(
          "dist" = "norm_mixture",
          "loglikelihood" = mixfit$loglik,
          "aic" = -2 * mixfit$loglik + 2 * mixfit.params,
          "bic" = -2 * mixfit$loglik + mixfit.params * log(length(x)),
          "mse" = mean(fit.residuals),
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
    if(fit.norm_mixture && !is.null(mixfit)) {
      mod$norm_mixture <- mixfit
      mod$norm_mixture$sd_scale <- sd.scale
    }
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
    x.adj <- (x - seg.start) + 1e-10
    # Extracting Model
    fti = models[[i]]
    # The bin size is now 1/sd_scale
    x.scale <- x.adj / fti$sd_scale

    if(class(fti) == "mixEM") {
      fit.data <- fti$x
      scalefactor <- length(fit.data) / fti$sd_scale
      distname <- "norm_mixture"
      dens <- dnorm_mixture(x.scale, fti)
    } else {
      fit.data <- fti$data
      # Adjust the scale of the segment
      # Multiple scale factor by this new bin size
      scalefactor <- length(fit.data) / fti$sd_scale

      # Generating data
      params <- c(as.list(fti$estimate), as.list(fti$fix.arg))
      distname <- fti$distname
      ddistname <- paste0("d", distname)
      call.params <- c(list(x = x.scale), as.list(params))
      dens <- do.call(ddistname, call.params)
    }
    dens.scale <- dens * scalefactor
    data.frame(
      "x" = x,
      "dens" = dens.scale,
      "col" = distname,
      row.names = NULL)
  }

  # Making a Nice Figure
  if(GENE %in% PARAMETERS$PLOT.MERGED.PEAKS) {

    distr.plotting.data = lapply(1:length(SEG.GR), comput.fti)
    # ggplot.plot(PARAMETERS=PARAMETERS, distr.plotting.data = distr.plotting.data, geneinfo=GENEINFO, bin.counts=BIN.COUNTS, seg.gr=SEG.GR, p=p)
    bpg.plot(PARAMETERS=PARAMETERS, distr.plotting.data = distr.plotting.data, geneinfo=GENEINFO, bin.counts=BIN.COUNTS, seg.gr=SEG.GR, p=p)

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
