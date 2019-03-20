library(VariantAnnotation)
library(tidyverse)

# path = "data/vcf/merged_exome.vcf.gz"
# path = "data/vcf/test_gt.vcf.gz"
# mindp = 4
read_gt <- function(path, mindp = 0, var_only = FALSE, tv_only = FALSE) {
  vcf <- readVcf(path)

  keep_pos <- elementNROWS(granges(vcf)$ALT) == 1

  if (var_only)
    keep_pos <- keep_pos & unlist(granges(vcf)$ALT) != ""

  gr_vcf <- granges(vcf)[keep_pos]

  info_df <- tibble(
      chrom = as.character(seqnames(gr_vcf)),
      pos = start(gr_vcf),
      REF = as.character(gr_vcf$REF),
      ALT = as.character(unlist(gr_vcf$ALT))
  )

  dp_mask <- geno(vcf)$DP[keep_pos, , drop = FALSE] %>% apply(2, function(i) ifelse(i >= mindp, i, NA))
  if ("chimp" %in% samples(header(vcf))) dp_mask[, "chimp"] <- 1

  gt_mat <- geno(vcf)$GT[keep_pos, , drop = FALSE] %>% replace(. == ".", NA) %>% replace(is.na(dp_mask), NA)
  mode(gt_mat) <- "numeric"
  gt_df <- gt_mat %>% as_tibble %>% mutate(reference = 0)

  # sanitize the sample names
  colnames(gt_df) <- str_replace_all(colnames(gt_df), "-", "_")

  res_df <- bind_cols(info_df, gt_df)

  if (tv_only) {
    res_df <- filter(res_df,  !((REF == "C" & ALT == "T") | (REF == "G" & ALT == "A")))
  }

  res_df
}
