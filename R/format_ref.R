#' Format the reference data.
#'
#' @param file_ref name of the reference data file
#' @param sumstats summary statistics data
#' @param group LD block group
#'
#' @return z-score and reference data
#' @export
#'
format_ref <- function(file_ref, sumstats, group = NULL) {
  message("Begin reading in reference data...")
  ref <- genio::read_plink(file_ref)
  message("The reference data set has ", nrow(ref$X), " lines.")

  if (!is.null(group)) {
    sumstats <- inner_join(sumstats, group, by = "SNP")
  }

  # merge the reference data
  ref_X <- bind_cols(SNP = rownames(ref$X), as_tibble(ref$X)) %>%
    drop_na() %>%
    semi_join(sumstats, by = "SNP")
  message("Remove missing data and merge the reference data with summary statistics data..., remaining ", nrow(ref_X), " SNPs.")

  # merge the summary statistics data with reference data
  sumstats <- semi_join(sumstats, ref_X, by = "SNP")
  message("Merge the summary statistics data with reference data..., remaining ", nrow(sumstats), " SNPs.")

  if (!is.null(group)) {
    z <- select(sumstats, Z, N, GROUP)
    X <- t(as.matrix(select(ref_X, -SNP)))
  } else {
    z <- select(sumstats, Z, N)
    X <- t(as.matrix(select(ref_X, -SNP)))
  }

  return(xpass_data = list(z = z, X = X))
}
