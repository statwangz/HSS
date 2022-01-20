#' Fit LD score regression.
#'
#' @param sumstats summary statistics data
#' @param ldsc LD score data
#' @param two_step logical
#' @param fix_intercept logical
#' @param n_blocks numeric
#' @param jackknife logical, whether to estimate the standard error
#'
#' @return results of LD score regression
#' @export
#'
ldsc_fit <- function(sumstats, ldsc,
                     two_step = T,
                     fix_intercept = F,
                     n_blocks = 200,
                     jackknife = T) {
  M <- nrow(ldsc) # number of SNPs included in the LD score estimation

  merged_data <- inner_join(sumstats, ldsc, by = "SNP")
  merged_data <- merged_data[with(merged_data, order(CHR, BP)), ]
  n_snps <- nrow(merged_data)

  if (two_step) {
    idx_step1 <- which(merged_data$CHI2 < 30)
  } else {
    idx_step1 <- NULL
  }

  # estimate heritability
  message("Begin LD score regression...")
  h2_se <- est_para(
    dat = select(merged_data, SNP, CHI2, N, L2),
    M = M,
    two_step = two_step,
    idx_step1 = idx_step1,
    fix_intercept = fix_intercept,
    n_blocks = n_blocks,
    jackknife = jackknife
  )

  return(h2_se)
}
