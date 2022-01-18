#' Estimate the parameters of LD score regression.
#'
#' @param dat data
#' @param M number of SNPs in LD score regression
#' @param two_step logical
#' @param idx_step1 vector, index of SNPs used in step 1
#' @param fix_intercept logical
#' @param n_blocks numeric
#'
#' @return slope and intercept in LD score regression
#' @export
#'
est_para <- function(dat, M,
                     two_step = T,
                     idx_step1 = NULL,
                     fix_intercept = F,
                     n_blocks = 200) {
  names(dat) <- c("SNP", "CHI2", "N", "L2")

  CHI2 <- dat$CHI2
  N <- dat$N
  L2 <- dat$L2
  L2 <- sapply(L2, function(x) {
    max(x, 1)
  })
  n_snps <- nrow(dat)

  # initial weights
  intercept <- 1
  h2_ini <- M * (mean(CHI2) - 1) / mean(L2 * N) # initial est of h2
  weights <- get_weights(h2_ini, intercept = 1, L2, N, M, L2)
  N_bar <- mean(N)
  x <- L2 * N / N_bar

  if (two_step) {
    # Two step estimator
    if (is.null(idx_step1)) {
      stop("please provide the index of SNPs in step 1!")
    }
    # step 1
    CHI2_s1 <- CHI2[idx_step1]
    L2_s1 <- L2[idx_step1]
    weights_s1 <- weights[idx_step1]
    N_s1 <- N[idx_step1]
    x_s1 <- x[idx_step1]
    seperator <- floor(seq(from = 1, to = length(idx_step1), length.out = (n_blocks + 1)))
    step1 <- irwls(CHi2_s1, L2_s1, update_x = x_s1, weights_s1, intercept = 1, M, N_s1, N_bar, fix_intercept, seperator)

    # step 2
    seperator_new <- idx_step1[seperator]
    step2 <- irwls(CHI2, L2, update_x = L2, weights, intercept = step1$intercept, M, N, N_bar, fix_intercept = T, new_seperator)
    h2 <- step2$h2
    intercept <- step1$intercept
  } else {
    seperator <- floor(seq(from = 1, to = length(CHI2), length.out = (n_blocks + 1)))
    step1 <- irwls(CHI2, L2, update_x = L2, weights, intercept = 1, M, N, N_bar, fix_intercept, seperator)
    h2 <- step1$h2
    intercept <- step1$intercept
  }

  return(list(h2 = h2, intercept = intercept))
}
