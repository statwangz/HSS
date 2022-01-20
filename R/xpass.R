#' Estimate heritability using XPASS
#'
#' @param z z-score
#' @param X reference genotypes matrix
#' @param K Kinship matrix
#' @param n sample size
#' @param Z covariates
#' @param group compute standard error by LD block Jackknife
#'
#' @return heritability and standard error
#' @export
#'
xpass <- function(z, X = NULL, K = NULL, n, Z = NULL, group = NULL) {
  if (is.null(X) & is.null(K)) {
    stop("please provide the information on reference genotypes data!")
  }

  if (is.null(K)) {
    X <- scale(X) / sqrt(ncol(X))
    K <- X %*% t(X)
  }

  m <- nrow(K)
  p <- length(z)

  # calculate h2 for y
  if (is.null(Z)) {
    Z <- matrix(1, m, 1)
    M <- diag(m) - matrix(1 / m, m, m)
    MK <- K
  } else {
    Z <- cbind(1, Z)
    M <- diag(m) - Z %*% solve(t(Z) %*% Z) %*% t(Z)
    MK <- M %*% K
  }

  q <- ncol(Z)

  trK <- sum(diag(MK))
  tr2K <- trK^2
  trK2 <- sum(MK^2)

  S <- (trK2 - tr2K / (m - q)) / (m - q)^2
  zz <- z^2 / n - 1 / n
  c <- sum(zz) / p

  h2 <- c / S

  if (is.null(group)) {
    # compute standard error by Jackknife
    c_jf <- (sum(zz) - zz) / (p - 1)
    se <- var(c_jf) / S / S * (p - 1)
  } else {
    group_num <- unique(group)
    ngroup <- length(group_num)
    zz_jf <- sapply(1:ngroup, function(j) {
      tmp <- sum(zz[group == group_num[j]])
      c(tmp, sum(group == group_num[j]))
    })
    c_jf <- (sum(zz) - zz_jf[1, ]) / (p - zz_jf[4, ])
    se <- var(c_jf) / S / S * (ngroup - 1)
  }

  return(list(h2 = h2, se = se))
}
