read_gt <- function(path) {
  vcf <- readVcf(path)

  gr_vcf <- granges(vcf)

  info_df <- tibble(
      chrom = levels(seqnames(gr_vcf)),
      pos = start(gr_vcf),
      REF = as.character(gr_vcf$REF),
      ALT = as.character(unlist(gr_vcf$ALT))
  )

  gt_mat <- geno(vcf)$GT %>% replace(. == ".", NA)
  mode(gt_mat) <- "numeric"
  gt_df <- as_tibble(gt_mat) %>% mutate(reference = 0)

  bind_cols(info_df, gt_df)
}
