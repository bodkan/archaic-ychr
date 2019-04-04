#' Read genotypes from a VCF file, returning a data frame object.
read_gt <- function(path, mindp = 0, var_only = FALSE, tv_only = FALSE, exclude = NA) {
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

  if (!is.na(exclude)) return(select(-one_of(exclude)))

  res_df
}


#' Read table of genotypes simulated by msprime.
read_simgt <- function(path) {
  suppressMessages(read_tsv(path)) %>%
    mutate(chrom = "simY", pos = round(pos), REF = "A", ALT = "T") %>%
    select(chrom, pos, REF, ALT, chimp = chimp0, everything())
}


#' Add a given proportion of sequencing/damage errors to a set of samples.
#' samples <- c("arch0", "arch1")
add_errors <- function(gt, rate, samples) {
  gt[, samples] <- map_dfr(gt[, samples], function(alleles) {
    mut_pos <- sort(sample(seq_along(alleles), size = length(alleles) * rate))
    alleles[mut_pos] <- as.integer(!as.logical(alleles[mut_pos]))
    alleles
  })
  gt
}
