format_ref <- function(file_ref, sumstats){
  ref <- genio::read_plink(file_ref)
  ref_snp <- select(ref$bim, id)
  names(ref_snp) <- c("SNP")
  
  # merge the summary statistics data
  sumstats <- inner_join(sumstats, ref_snp, by = "SNP")
  message("Merge the summary statistics with reference data..., remaining ", nrow(sumstats), " SNPs.")
  
  # merge the reference data
  ref_X <- bind_cols(SNP = rownames(ref$X), as_tibble(ref$X)) %>%
    inner_join(select(sumstats, SNP)) %>%
    select(-SNP)
  
  z <- select(sumstats, Z)
  
  return(xpass_data <- list(z, ref_X))
}