#' Group the SNPs using LD block data.
#'
#' @param ref_file name of the reference data file
#' @param snp_col name of a column with rs numbers. The default is `id`.
#' @param pos_col name of a column with position. The default is `pos`.
#'
#' @return group tibble
#' @export
#'
ld_block_group <- function(ref_file, snp_col = "id", pos_col = "pos") {
  bim <- genio::read_plink(ref_file)$bim
  pos <- select(bim, one_of(pos_col))
  snp <- select(bim, one_of(snp_col))

  data(EUR_ld_block)
  ngroup <- nrow(block)

  group <- apply(pos,
    MARGIN = 1, FUN = pos_group,
    start = block$start, stop = block$stop, ngroup = ngroup
  )

  group <- tibble(snp, group)
  names(group) <- c("SNP", "GROUP")

  return(group)
}
