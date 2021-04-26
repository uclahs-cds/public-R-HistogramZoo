
calculate.density = function(
m,
x = NULL,
seg.start,
seg.end,
stepsize = 1,
scale.density = T,
return.df = F
){

  # Range of x
  if(is.null(x)){
    x.range = seq(seg.start, seg.end, stepsize)
    x.range.adjusted <- (x.range - seg.start) + 1e-10
  } else {
    x.range.adjusted = x
  }

  # Calculate Density
  if(class(m) == "mixEM") {
    m.data <- m$x
    scalefactor <- length(m.data)
    distname <- "norm_mixture"
    dens <- dnorm_mixture(x.range.adjusted, m)
  } else {
    m.data <- m$data
    scalefactor <- length(m.data)
    params <- c(as.list(m$estimate), as.list(m$fix.arg))
    distname <- m$distname
    ddistname <- paste0("d", distname)
    call.params <- c(list(x = x.range.adjusted), as.list(params))
    dens <- do.call(ddistname, call.params)
  }

  # Scaling Density
  if(scale.density){
    dens <- dens * scalefactor
  }

  # Return format
  if(return.df){
    data.frame(
      "x" = x.range,
      "dens" = dens,
      "dist" = distname,
      row.names = NULL)
  } else {
    return(dens)
  }
}
