#' Estimate the parameters of LD score regression.
#'
#' @param dat data
#' @param M number of SNPs in LD score regression
#' @param two_step logical
#' @param idx_step1 vector, index of SNPs used in step 1
#' @param fix_intercept logical
#' @param n_blocks numeric
#' @param jackknife logical, whether to estimate the standard error
#'
#' @return slope and intercept in LD score regression
#' @export
#'
est_para <- function(dat, M,
                     two_step = T,
                     idx_step1 = NULL,
                     fix_intercept = F,
                     n_blocks = 200,
                     jackknife = T) {
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

  if (fix_intercept) {
    two_step <- F
  }
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
    step1 <- irwls(CHI2_s1, L2_s1, update_x = x_s1, weights_s1, intercept = 1, M, N_s1, N_bar, fix_intercept, seperator, jackknife = jackknife)

    # step 2
    seperator_new <- idx_step1[seperator]
    step2 <- irwls(CHI2, L2, update_x = L2, weights, intercept = step1$intercept, M, N, N_bar, fix_intercept = T, seperator_new, jackknife = jackknife)
    h2 <- step2$h2
    intercept <- step1$intercept

    if (jackknife) {
      ## cpmbine step 1 and step 2
      c <- sum(weights^2 * x) / sum(weights^2 * x^2)

      est <- c(step2$jk_est, step1$intercept)
      delete_values <- matrix(0, nrow = n_blocks, ncol = 2)
      delete_values[, 2] <- step1$delete_values[, 2]
      delete_values[, 1] <- (step2$delete_values - c * (step1$delete_values[, 2] - step1$intercept))

      pseudo_values <- matrix(n_blocks * est, nrow = n_blocks, ncol = 2, byrow = T) - (n_blocks - 1) * delete_values

      jackknife_cov <- cov(pseudo_values) / n_blocks
      jackknife_se <- sqrt(diag(jackknife_cov))

      coef_cov <- jackknife_cov[1, 1] / (N_bar^2) * M^2
      h2_se <- sqrt(coef_cov)
      intercept_se <- jackknife_se[2]
    }
  } else {
    seperator <- floor(seq(from = 1, to = length(CHI2), length.out = (n_blocks + 1)))
    step1 <- irwls(CHI2, L2, update_x = L2, weights, intercept = 1, M, N, N_bar, fix_intercept, seperator, jackknife = jackknife)
    h2 <- step1$h2
    intercept <- step1$intercept

    if (jackknife) {
      h2_se <- step1$h2_se
      intercept_se <- step1$intercept_se
      delete_values <- step1$delete_values
    }
  }

  if (jackknife) {
    return(list(
      h2 = h2, h2_se = h2_se,
      intercept = intercept,
      intercept_se = intercept_se,
      delete_values = delete_values[, 1]
    ))
  } else {
    return(list(h2 = h2, intercept = intercept))
  }
}
