# Peak method from https://github.com/stas-g/findPeaks/blob/master/find_peaks.R
# https://stats.stackexchange.com/questions/22974/how-to-find-local-peaks-valleys-in-a-series-of-data/164830#164830
#
# plot(BIN.COUNTS$start, BIN.COUNTS$Coverage, type = "s")
# p <- find.peaks(-BIN.COUNTS$Coverage, m = 150)
# points(BIN.COUNTS$start[p], BIN.COUNTS$Coverage[p], col = 'red')
find.peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}
