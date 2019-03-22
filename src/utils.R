library(VariantAnnotation)
library(tidyverse)

# path = "data/vcf/merged_exome.vcf.gz"
# path = "data/vcf/test_gt.vcf.gz"
# mindp = 4
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



#
# Generate table of sample names and their population assignments and ages
#
sample_info <- function(gt) {
    emh <- tibble(
        name = c("ustishim", "bichon", "kk1", "loschbour", "mota"),
        age = c(45000, 13665, 9720, 8050, 4500),
        pop = "EMH"
    )

    modern <- tibble(
        name = c(str_subset(colnames(gt), "^S_"), "reference", "a00"), #, "a00_1", "a00_2"),
        pop = case_when(
            name == "reference" ~ "reference",
            name %in% c("S_Burmese_1", "S_Thai_1", "S_Han_2", "S_Dai_2", "S_Punjabi_1", "S_Papuan_2", "S_Karitiana_1") ~ "EastEur",
            name %in% c("S_BedouinB_1", "S_Turkish_1", "S_French_1", "S_Finnish_2", "S_Sardinian_1", "S_Saami_2") ~ "WestEur",
            name %in% c("ustishim", "bichon", "kk1", "loschbour", "mota") ~ "EMH",
            TRUE ~ "Africa"
        ),
        age = 0
    )

    bind_rows(emh, modern)
}



#
# Assign mutations to all possible branches of a tree with 4 leaves.
#
site_patterns <- function(df, w, x, y, z) {
    tibble(
        a = as.integer(df[[w]] == df[[x]] & df[[y]] == df[[z]] & df[[w]] != df[[y]]),
        b = as.integer(df[[w]] == df[[z]] & df[[x]] == df[[y]] & df[[w]] != df[[x]]),
        c = as.integer(df[[w]] == df[[y]] & df[[x]] == df[[z]] & df[[w]] != df[[z]]),
        d = as.integer(df[[w]] == df[[x]] & df[[w]] == df[[z]] & df[[w]] != df[[y]]),
        e = as.integer(df[[w]] == df[[x]] & df[[w]] == df[[y]] & df[[w]] != df[[z]]),
        f = as.integer(df[[w]] == df[[y]] & df[[w]] == df[[z]] & df[[w]] != df[[x]])
    )
}



#
# Count the occurences of different site patterns.
#
sum_patterns <- function(gt, w, x, y, z) {
    df <- gt[, c(w, x, y, z)] %>% .[complete.cases(.), ]
    site_patterns(df, w, x, y, z) %>% summarise_all(sum) %>% mutate(total = nrow(df))
}




