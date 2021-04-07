# Peak method from https://github.com/stas-g/findPeaks/blob/master/find_peaks.R
# https://stats.stackexchange.com/questions/22974/how-to-find-local-peaks-valleys-in-a-series-of-data/164830#164830
#
# plot(BIN.COUNTS$start, BIN.COUNTS$Coverage, type = "s")
# p <- find.peaks(-BIN.COUNTS$Coverage, m = 150)
# points(BIN.COUNTS$start[p], BIN.COUNTS$Coverage[p], col = 'red')
find.peaks <- function(x, m = 3, diff.threshold = 0.001){
  diff.x = diff(x, na.pad = FALSE)
  diff.x[abs(diff.x) < diff.threshold] <- 0
  shape <- diff(sign(diff.x))
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


# Moving Average From https://stackoverflow.com/questions/743812/calculating-moving-average
moving.average = function(x, n = 5) {
  stats::filter(x, rep(1 / n, n), sides = 2)
}

# Adapted from exomepeak
remove.local.abnormalities <- function(bin.counts, max_background_fold_increase=2, background_window=10){

  # get position
  pos_table = bin.counts$Coverage
  pos_mapped = bin.counts$start
  no_pos_mapped=length(pos_mapped)

  # prepare new table
  new_table=bin.counts$Coverage

  # filter
  for (i in 1:no_pos_mapped) {
    ID=which(abs(pos_mapped[i]-pos_mapped) < (background_window/2))
    else_count = sum(pos_table[ID])-sum(pos_table[i])
    max_background = round(else_count*max_background_fold_increase/background_window)
    new_table[i] = max(1,min(pos_table[i],max_background))
  }

  return(new_table)
}
