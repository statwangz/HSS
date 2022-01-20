pos_group <- function(pos, start, stop, ngroup) {
  group <- 0
  i <- 1

  while (i <= ngroup) {
    if (pos >= start[i] & pos <= stop[i]) {
      group <- i
      break
    }
    i <- i + 1
  }

  return(group)
}
