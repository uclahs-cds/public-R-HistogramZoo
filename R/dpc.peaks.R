.dpc.peaks = function(GENEPEAKSGR, PARAMETERS){

  # Sampling from the set of peaks
  TILED.PEAKS.GR = unlist(GenomicRanges::tile(GENEPEAKSGR, width = PARAMETERS$DP.RESOLUTION))

  # Normalizing Data
  startvec = GenomicRanges::start(TILED.PEAKS.GR)
  startvec.mean = mean(startvec)
  startvec.sd = sd(startvec)
  # This is for when there's 1 short peak or all samples have the exact same peak (Fix this at some point)
  startvec.sd = ifelse( is.na(startvec.sd) | startvec.sd == 0, 1, startvec.sd)
  startvec.scaled = (startvec - startvec.mean)/startvec.sd

  # Running Dirichlet Process
  set.seed(PARAMETERS$DP.SEED)
  dp = dirichletprocess::DirichletProcessGaussian(
    y = startvec.scaled,
    g0Priors = c(0, 1, 1, 1),
    alphaPriors = PARAMETERS$DP.ALPHA.PRIORS)

  dp = dirichletprocess::Fit(
    dpObj = dp,
    its = PARAMETERS$DP.ITERATIONS,
    updatePrior = F,
    progressBar = F)

  # Generating a data.frame from the dp object
  dp_data = data.frame(
    "Weights" = dp$weights,
    "Mu" = c(dp$clusterParameters[[1]]),
    "Sigma" = c(dp$clusterParameters[[2]]),
    stringsAsFactors = F)
  dp_data$Mu = (dp_data$Mu*startvec.sd)+startvec.mean
  dp_data$Sigma = (dp_data$Sigma*startvec.sd)

  return(list("dp" = dp, "dp_data" = dp_data, "startvec" = startvec, "startvec.scaled" = startvec.scaled, "startvec.mean" = startvec.mean, "startvec.sd" = startvec.sd))

}
