library(VariantAnnotation)
library(tidyverse)

# path = "data/vcf/merged_full.vcf.gz"
# mindp = 4
read_gt <- function(path, mindp = 0) {
  vcf <- readVcf(path)

  gr_vcf <- granges(vcf) %>% .[elementNROWS(.$ALT) == 1]

  info_df <- tibble(
      chrom = as.character(seqnames(gr_vcf)),
      pos = start(gr_vcf),
      REF = as.character(gr_vcf$REF),
      ALT = as.character(unlist(gr_vcf$ALT))
  )

  dp_mask <- geno(vcf)$DP %>% apply(2, function(i) ifelse(i >= mindp, i, NA))
  if ("chimp" %in% samples(header(vcf))) dp_mask[, "chimp"] <- 1

  gt_mat <- geno(vcf)$GT %>% replace(. == ".", NA) %>% replace(is.na(dp_mask), NA)
  mode(gt_mat) <- "numeric"
  gt_df <- gt_mat %>% as_tibble %>% mutate(reference = 0)

  bind_cols(info_df, gt_df)
}
