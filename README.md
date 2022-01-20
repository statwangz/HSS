
<!-- README.md is generated from README.Rmd. Please edit that file -->

# HSS

<!-- badges: start -->
<!-- badges: end -->

The R package `HSS` implements two methods `LD score regression` and
`XPASS` to estimate heritability only using summary statistics data.

## Installation

You can install the development version of HSS from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("statwangz/HSS")
```

## Example

Here is a basic example:

``` r
# Download data
# Download BMI summary statistics from https://atlas.ctglab.nl/
# Download HapMap 3 data from Broad Institute
# wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
# Download LD score data from Broad Institute
# wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
# Download 1000 Genomes data from https://www.cog-genomics.org/plink/2.0/resources

library(HSS)

file_sumstats <- "f.23104.0.0_res.EUR.sumstats.MACfilt.txt"
# read in HapMap 3 data
w_hm3.snplist <- readr::read_delim("w_hm3.snplist.bz2",  "\t",
                                   escape_double = F, trim_ws = T, progress = T)
# load MHC region SNPs, CHR = 6 & BP > 28000000 & BP < 34000000
data(snps_mhc)
# format BMI summary statistics data
BMI_sumstats <- format_sumstats(file_sumstats,
                                snps_merge = w_hm3.snplist,
                                snps_mhc = snps_mhc,
                                snp_col = "SNPID_UKB",
                                beta_col = "BETA",
                                se_col = "SE",
                                freq_col = "MAF_UKB",
                                a1_col = "A1",
                                a2_col = "A2",
                                p_col = "P",
                                n_col = "NMISS",
                                info_col = "INFO_UKB")
head(BMI_sumstats, n = 5)

# format LD score data
file_ldsc <- "eur_w_ld_chr"
ldsc <- format_ldsc(file_ldsc)
head(ldsc, n = 5)

# LD score regression
res_ldsc <- ldsc_fit(BMI_sumstats, ldsc)
cat("The estimate of heritability is ", round(res_ldsc$h2, 2), ".\n", sep = "")

# format reference data
file_ref <- "1000G.EUR.QC.hm3.ind"
xpass_data <- format_ref(file_ref, BMI_sumstats)
z <- xpass_data$z
X <- xpass_data$X

# XPASS
h2 <- xpass(z = pull(z, Z), K = K, n = pull(z, N))
cat("The estimate of heritability is ", round(h2, 2), ".\n", sep = "")
```

Please see
[here](https://github.com/statwangz/HSS/tree/main/ten_traits_analysis)
for more examples which analyze other traits.

## Development

The `HSS` package is developed by Zhiwei Wang
(<zhiwei.wang@connect.ust.hk>).

## Contact information

Please feel free to contact Zhiwei Wang (<zhiwei.wang@connect.ust.hk>)
if any enquiry.

## References

Bulik-Sullivan, B.K., Loh, P.-R., Finucane, H.K., Ripke, S., Yang, J.,
Patterson, N., Daly, M.J., Price, A.L., Neale, B.M., 2015. LD Score
regression distinguishes confounding from polygenicity in genome-wide
association studies. Nat Genet 47, 291–295.
<https://doi.org/10.1038/ng.3211>

Cai, M., Xiao, J., Zhang, S., Wan, X., Zhao, H., Chen, G., Yang, C.,
2021. A unified framework for cross-population trait prediction by
leveraging the genetic correlation of polygenic traits. The American
Journal of Human Genetics 108, 632–655.
<https://doi.org/10.1016/j.ajhg.2021.03.002>

Hu, X., Zhao, J., Lin, Z., Wang, Y., Peng, H., Zhao, H., Wan, X., Yang,
C., 2021. MR-APSS: a unified approach to Mendelian Randomization
accounting for pleiotropy and sample structure using genome-wide summary
statistics. <https://doi.org/10.1101/2021.03.11.434915>
