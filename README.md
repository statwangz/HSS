
<!-- README.md is generated from README.Rmd. Please edit that file -->

# HSS

<!-- badges: start -->
<!-- badges: end -->

The goal of HSS is to estimate heritability only using summary
statistics data.

## Installation

You can install the development version of HSS from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("statwangz/HSS")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
# Download data

# Download BMI summary statistics from https://atlas.ctglab.nl/

# Download HapMap 3 data from Broad Institute
# wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2

# Download LD score data from Broad Institute
# wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2

library(HSS)
#> Loading required package: tidyverse
#> ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──
#> ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
#> ✓ tibble  3.1.6     ✓ dplyr   1.0.7
#> ✓ tidyr   1.1.4     ✓ stringr 1.4.0
#> ✓ readr   2.1.1     ✓ forcats 0.5.1
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> x dplyr::filter() masks stats::filter()
#> x dplyr::lag()    masks stats::lag()
```

## Development

The HSS package is developed by Zhiwei Wang
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

Hu, X., Zhao, J., Lin, Z., Wang, Y., Peng, H., Zhao, H., Wan, X., Yang,
C., 2021. MR-APSS: a unified approach to Mendelian Randomization
accounting for pleiotropy and sample structure using genome-wide summary
statistics. <https://doi.org/10.1101/2021.03.11.434915>

Cai, M., Xiao, J., Zhang, S., Wan, X., Zhao, H., Chen, G., Yang, C.,
2021. A unified framework for cross-population trait prediction by
leveraging the genetic correlation of polygenic traits. The American
Journal of Human Genetics 108, 632–655.
<https://doi.org/10.1016/j.ajhg.2021.03.002>
