
find_new_segments = function(gaps.df) {
  new.segments = c()
  for(i in rownames(gaps.df)) {
    r = gaps.df[i, ]
    # Have a large gap
    if(r$end - r$start > 1) {
      left.diff = abs(r$seg.start - r$start)
      right.diff = abs(r$seg.end-r$end)
      if(left.diff < right.diff && right.diff > 2) {
        new.segments = c(new.segments, r$end - 1)
      } else if(left.diff > right.diff && left.diff > 2){
        new.segments = c(new.segments, r$start + 1)
      }
    }
  }
  return(new.segments)
}
