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



read_simgt <- function(path) {
  suppressMessages(read_tsv(path)) %>%
    mutate(chrom = "simY", pos = round(pos), REF = "A", ALT = "T") %>%
    select(chrom, pos, REF, ALT, chimp = chimp0, everything())
}



read_siminfo <- function(gt) {
  emh <- tibble(name = str_subset(colnames(gt), "ustishim"), age = 45000, pop = "EMH")
  afr <- tibble(name = str_subset(colnames(gt), "afr"), age = 0, pop = "Africa")
  westeur <- tibble(name = str_subset(colnames(gt), "eur"), age = 0, pop = "WestEur")
  easteur <- tibble(name = str_subset(colnames(gt), "asn"), age = 0, pop = "EastEur")
  bind_rows(emh, afr, westeur, easteur)
}



#
# Generate table of sample names and their population assignments and ages
#
read_info <- function(gt) {
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
    n <- if (gt$chrom[1] == "simY") max(gt$pos) else nrow(df)
    site_patterns(df, w, x, y, z) %>% summarise_all(sum) %>% mutate(total = n)
}



# Calculate T_MRCA of all pairs of Africans and non-Africans, using a set
# of EMHs of known age for branch shortening-based mutation rate estimation.
calculate_tAfrRef <- function(gt, samples) {
    refs <- filter(samples, pop != "Africa", pop != "EMH")$name
    afrs <- filter(samples, pop == "Africa")$name
    emhs <- filter(samples, pop == "EMH")$name

    map_dfr(refs, function(ref) {
        map_dfr(afrs, function(afr) {
            map_dfr(emhs, function(emh) {
                sum_patterns(gt, w = "chimp", x = afr, y = ref, z = emh) %>%
                    mutate(ref = ref, afr = afr, emh = emh)
            })
        })
    }) %>%
    inner_join(samples, by = c("emh" = "name")) %>%
    mutate(
        muts_per_year = (d - e) / age,
        mut_rate = muts_per_year / total,
        tmrca_ad = (a + d) / (mut_rate * total),
        tmrca_f = f / (mut_rate * total),
        tmrca_avg = (tmrca_ad + tmrca_f) / 2
    )
}





# Calculate T_MRCA of each of the specified archaic individuals
# with modern human Y chromosomes.
calculate_tArchRef <- function(gt, samples, tAfrRef) {
    refs <- filter(samples, pop != "Africa", pop != "EMH")$name
    afrs <- filter(samples, pop == "Africa")$name
    archaics <- colnames(select(gt, -c(chrom, pos, REF, ALT, chimp), -one_of(samples$name)))

    tAfrRef_ui <- filter(tAfrRef, emh == "ustishim") %>%
        select(afr, ref, mut_rate, tmrca_ad, tmrca_f, tmrca_avg)

    map_dfr(refs, function(ref) {
        map_dfr(afrs, function(afr) {
            map_dfr(archaics, function(arch) {
                sum_patterns(gt, w = "chimp", x = arch, y = ref, z = afr) %>%
                    mutate(arch = arch, afr = afr, ref = ref)
            })
        })
    }) %>% inner_join(tAfrRef_ui, by = c("afr", "ref")) %>%
    mutate(p = a / (a + d + e),
           alpha = (1 + p) / (1 - p),
           tmrca_arch = tmrca_avg * alpha) %>%
    select(arch, afr, ref, tmrca_arch, alpha, mut_rate, tmrca_ad, tmrca_f, tmrca_avg, everything())
}
