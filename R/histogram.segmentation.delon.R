
# For the Grenador estimator, Pool Adjacent Violators algorithm
library(isotone)

# Functions ---------------------------------------------------------------

# Observed distribution
rab = function(N, obs, a, b, seg.size){
  
  obs = obs[obs >= a & obs < b]
  bk.pts = seq(a, b, seg.size)
  dens = table(cut(obs, breaks = bk.pts))/N
  dens = as.numeric(dens)
  return(dens)
}

# Grenander estimator, Pool Adjacent Violators
pab = function(dens, increasing = T){
  
  if(increasing){
    est = isotone::gpava(z = 1:length(dens), y = dens)
    est = est$x
  } else {
    est = isotone::gpava(z = 1:length(dens), y = rev(dens))
    est = rev(est$x)
  }
  
  return(est)
}

# NFA
nfa = function(rab., pab.){
  
  L = length(rab.)
  coeff = L*(L+1)/2
  
  bnkp = c()
  for(i in 1:length(rab.)){
    rab.i = rab.[i]
    pab.i = pab.[i]
    
    if(rab.i >= pab.i){
      bnkp[i] = pbinom(N*rab.i, N, pab.i, lower.tail = F)
    } else {
      bnkp[i] = pbinom(N*(1-rab.i), N, (1-pab.i), lower.tail = F)
    }
  }
  bnkp = prod(bnkp)*coeff
  
  return(bnkp)
}

# Testing for satisfaction of the unimodal hypothesis
unimodal.hypothesis = function(rab. , epsilon = 1){
  
  # Initializing
  unimodal.binary = F
  max.point = order(-rab.)
  i = 1
  thres = epsilon/2

    while(!isTRUE(unimodal.binary) | i < length(max.point)){
    
    # Choosing midpoint to test
    midpoint.to.test = max.point[i]
    pab.increasing = pab(rab.[1:midpoint.to.test], increasing = T)
    pab.decreasing = pab(rab.[midpoint.to.test:length(rab)], increasing = F)
    
    # NFA
    nfa.increasing = nfa(rab. = rab.[1:midpoint.to.test], pab. = pab.increasing)
    nfa.decreasing = nfa(rab. = rab.[midpoint.to.test:length(rab)], pab. = pab.decreasing)
    
    # Updating
    unimodal.binary = (nfa.increasing < thres & nfa.decreasing < thres)
    i = i+1
    } 
  
  return(unimodal.binary)
  
}

# Fine to Coarse Segmentation Algorithm
ftc.segments = function(obs, seg.size){
  
  # Initializing segmentation
  a = floor(min(obs))
  b = ceiling(max(obs))
  bins = seq(a, b, seg.size)
  
  # Generating density vector
  rab. = rab(N = length(obs), obs = obs, a = a, b = b, seg.size = seg.size)
  
  # counter, this checks whether every test has been done
  ct = 0
  
  while(ct < (length(bins)-2)){
    
    # Sampling a point
    i = sample(2:(length(bins)-1), 1)
    test.i = unimodal.hypothesis(rab.[bins[i-1]:bins[i+1]])
    
    # Updating bins
    if(test.i){
      bins = bins[!bins %in% bins[i]]
      ct = 0
    } else {
      ct = ct+1
    }
  }
  
}


# Testing -----------------------------------------------------------------

set.seed(314)
simulated.data = c(rnorm(100, 5, 1), rnorm(100, 12, 3))
hist(simulated.data, breaks = 1:20)
summary(simulated.data)

# Make a plot

