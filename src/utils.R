library(VariantAnnotation)
library(tidyverse)

# path = "data/vcf/merged_full.vcf.gz"
# path = "data/vcf/test_gt.vcf.gz"
# mindp = 4
read_gt <- function(path, mindp = 0, var_only = FALSE) {
  vcf <- readVcf(path)

  bialellic_pos <- elementNROWS(granges(vcf)$ALT) == 1

  if (var_only)
    bialellic_pos <- bialellic_pos & unlist(granges(vcf)$ALT) != ""

  gr_vcf <- granges(vcf)[bialellic_pos]

  info_df <- tibble(
      chrom = as.character(seqnames(gr_vcf)),
      pos = start(gr_vcf),
      REF = as.character(gr_vcf$REF),
      ALT = as.character(unlist(gr_vcf$ALT))
  )

  dp_mask <- geno(vcf)$DP %>% .[bialellic_pos, ] %>% apply(2, function(i) ifelse(i >= mindp, i, NA))
  if ("chimp" %in% samples(header(vcf))) dp_mask[, "chimp"] <- 1

  gt_mat <- geno(vcf)$GT %>% .[bialellic_pos, ] %>% replace(. == ".", NA) %>% replace(is.na(dp_mask), NA)
  mode(gt_mat) <- "numeric"
  gt_df <- gt_mat %>% as_tibble %>% mutate(reference = 0)

  bind_cols(info_df, gt_df)
}
