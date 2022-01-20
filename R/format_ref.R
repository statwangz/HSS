#' Format the reference data.
#'
#' @param X_ref reference genotypes data
#' @param sumstats summary statistics data
#' @param group LD block group
#'
#' @return z-score and reference data
#' @export
#'
format_ref <- function(X_ref, sumstats, group = NULL) {
  message("The reference data set has ", nrow(X_ref), " lines.")

  if (!is.null(group)) {
    message("Group the SNPs in summary statistics data...")
    sumstats <- inner_join(sumstats, group, by = "SNP")
  }

  # merge the reference data
  X_ref <- bind_cols(SNP = rownames(X_ref), as_tibble(X_ref)) %>%
    drop_na() %>%
    semi_join(sumstats, by = "SNP")
  message("Remove missing data and merge the reference data with summary statistics data..., remaining ", nrow(X_ref), " SNPs.")

  # merge the summary statistics data with reference data
  sumstats <- semi_join(sumstats, X_ref, by = "SNP")
  message("Merge the summary statistics data with reference data..., remaining ", nrow(sumstats), " SNPs.")

  if (!is.null(group)) {
    z <- select(sumstats, Z, N, GROUP)
    X <- t(as.matrix(select(X_ref, -SNP)))
  } else {
    z <- select(sumstats, Z, N)
    X <- t(as.matrix(select(X_ref, -SNP)))
  }

  return(xpass_data = list(z = z, X = X))
}
