data(trait_info)
trait_name <- trait_info$Trait
trait_name

file_name <- paste0(trait_name, ".EUR.sumstats.txt.gz")
format_name <- paste0("sumstats/", trait_name, "_sumstats.txt")

# read in HapMap 3 data
w_hm3.snplist <- readr::read_delim("w_hm3.snplist.bz2",  "\t",
                                   escape_double = F, trim_ws = T, progress = T)
# load MHC region SNPs, CHR = 6 & BP > 28000000 & BP < 34000000
data(snps_mhc)

Angina_sumstats <- format_sumstats(file_name[1],
                                   snps_merge = w_hm3.snplist,
                                   snps_mhc = snps_mhc,
                                   snp_col = "SNPID_UKB",
                                   or_col = "OR",
                                   se_col = "SE",
                                   freq_col = "MAF_UKB",
                                   a1_col = "A1",
                                   a2_col = "A2",
                                   p_col = "P",
                                   n_col = "NMISS",
                                   info_col = "INFO_UKB")
write_delim(Angina_sumstats, format_name[1])
rm(Angina_sumstats)

BMI_sumstats <- format_sumstats(file_name[2],
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
write_delim(BMI_sumstats, format_name[2])
rm(BMI_sumstats)

Depression_sumstats <- format_sumstats(file_name[3],
                                       snps_merge = w_hm3.snplist,
                                       snps_mhc = snps_mhc,
                                       snp_col = "RSID",
                                       freq_col = "MAF_UKB",
                                       a1_col = "A1",
                                       a2_col = "A2",
                                       z_col = "Z",
                                       p_col = "P",
                                       n_col = "N",
                                       info_col = "INFO_UKB")
write_delim(Depression_sumstats, format_name[3])
rm(Depression_sumstats)

Hair_sumstats <- format_sumstats(file_name[4],
                                 snps_merge = w_hm3.snplist,
                                 snps_mhc = snps_mhc,
                                 snp_col = "SNPID_UKB",
                                 or_col = "OR",
                                 se_col = "SE",
                                 freq_col = "MAF_UKB",
                                 a1_col = "A1",
                                 a2_col = "A2",
                                 p_col = "P",
                                 n_col = "NMISS",
                                 info_col = "INFO_UKB")
write_delim(Hair_sumstats, format_name[4])
rm(Hair_sumstats)

HBP_sumstats <- format_sumstats(file_name[5],
                                snps_merge = w_hm3.snplist,
                                snps_mhc = snps_mhc,
                                snp_col = "SNPID_UKB",
                                or_col = "OR",
                                se_col = "SE",
                                freq_col = "MAF_UKB",
                                a1_col = "A1",
                                a2_col = "A2",
                                p_col = "P",
                                n_col = "NMISS",
                                info_col = "INFO_UKB")
write_delim(HBP_sumstats, format_name[5])
rm(HBP_sumstats)

Height_sumstats <- format_sumstats(file_name[6],
                                   snps_merge = w_hm3.snplist,
                                   snps_mhc = snps_mhc,
                                   snp_col = "SNPID_UKB",
                                   beta_col = "BETA",
                                   se_col = "SE",
                                   freq_col = "MAF",
                                   a1_col = "A1",
                                   a2_col = "A2",
                                   p_col = "P",
                                   n_col = "NMISS",
                                   info_col = "INFO_UKB")
write_delim(Height_sumstats, format_name[6])
rm(Height_sumstats)

Income_sumstats <- format_sumstats(file_name[7],
                                   snps_merge = w_hm3.snplist,
                                   snps_mhc = snps_mhc,
                                   snp_col = "SNPID_UKB",
                                   beta_col = "BETA",
                                   se_col = "SE",
                                   freq_col = "MAF",
                                   a1_col = "A1",
                                   a2_col = "A2",
                                   p_col = "P",
                                   n_col = "NMISS",
                                   info_col = "INFO_UKB")
write_delim(Income_sumstats, format_name[7])
rm(Income_sumstats)

Intelligence_sumstats <- format_sumstats(file_name[8],
                                         snps_merge = w_hm3.snplist,
                                         snps_mhc = snps_mhc,
                                         snp_col = "SNPID_UKB",
                                         beta_col = "BETA",
                                         se_col = "SE",
                                         freq_col = "MAF",
                                         a1_col = "A1",
                                         a2_col = "A2",
                                         p_col = "P",
                                         n_col = "NMISS",
                                         info_col = "INFO_UKB")
write_delim(Intelligence_sumstats, format_name[8])
rm(Intelligence_sumstats)

MDD_sumstats <- format_sumstats(file_name[9],
                                snps_merge = w_hm3.snplist,
                                snps_mhc = snps_mhc,
                                snp_col = "SNPID_UKB",
                                or_col = "OR",
                                se_col = "SE",
                                freq_col = "MAF",
                                a1_col = "A1",
                                a2_col = "A2",
                                p_col = "P",
                                n_col = "NMISS",
                                info_col = "INFO_UKB")
write_delim(MDD_sumstats, format_name[9])
rm(MDD_sumstats)

Tanning_sumstats <- format_sumstats(file_name[10],
                                    snps_merge = w_hm3.snplist,
                                    snps_mhc = snps_mhc,
                                    snp_col = "SNPID_UKB",
                                    or_col = "BETA",
                                    se_col = "SE",
                                    freq_col = "MAF",
                                    a1_col = "A1",
                                    a2_col = "A2",
                                    p_col = "P",
                                    n_col = "NMISS",
                                    info_col = "INFO_UKB")
write_delim(Tanning_sumstats, format_name[10])
rm(Tanning_sumstats)
