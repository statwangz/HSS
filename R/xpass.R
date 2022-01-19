#' Estimate heritability using XPASS
#'
#' @param z z-score
#' @param X reference genotypes matrix
#' @param K Kinship matrix
#' @param n sample size
#' @param Z covariates
#'
#' @return heritability
#' @export
#'
xpass <- function(z, X = NULL, K = NULL, n, Z = NULL) {
  if(is.null(X) & is.null(K)){
    stop("please provide the information on reference genotypes data!")
  }

  if(is.null(K)){
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

  return(h2)
}
