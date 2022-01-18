get_weights <- function(h2, intercept, L2, N, M, x) {
  h2 <- max(h2, 0)
  h2 <- min(h2, 1)
  wld <- sapply(L2, function(x) {
    max(x, 1)
  })
  ld <- sapply(x, function(x) {
    max(x, 1)
  })
  c <- N * h2 / M
  het_w <- 1 / (2 * (intercept + c * ld)^2)
  # cat("mean het_w ", mean(het_w),"\n")
  oc_w <- 1 / wld
  # cat("mean oc_w ", mean(oc_w), "\n")
  w <- sqrt(het_w * oc_w)
  return(w)
}
