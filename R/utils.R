#' Generate table of sample names and their population assignments and ages
#' from the simulated data.
read_siminfo <- function(gt) {
  emh <- tibble(name = str_subset(colnames(gt), "ustishim"), age = 45000, pop = "EMH")
  afr <- tibble(name = str_subset(colnames(gt), "afr"), age = 0, pop = "Africa")
  westeur <- tibble(name = str_subset(colnames(gt), "eur"), age = 0, pop = "WestEur")
  easteur <- tibble(name = str_subset(colnames(gt), "asn"), age = 0, pop = "EastEur")
  bind_rows(emh, afr, westeur, easteur)
}


#' Generate table of sample names and their population assignments and ages.
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


#' Calculate proportion of different classes of SNPs observed in each sample.
#' @param normalize Normalize counts by some SNP class? Either FALSE or "T-C", etc.
snp_props <- function(gt, normalize = NA, remove = NA) {
  snp_counts <- gt %>%
    mutate(snp = REF %+% "-" %+% ALT) %>%
    filter(!snp %in% remove) %>%
    select(-(chrom:ALT)) %>%
    group_by(snp) %>%
    summarise_all(~sum(., na.rm = TRUE))

  total_sites <-
    gt %>%
    select(-(chrom:ALT)) %>%
    summarise_all(~ sum(!is.na(.)))

  snp_props <-
    snp_counts[, -1] %>%
    colnames %>%
    map(~ snp_counts[[.x]] / total_sites[[.x]]) %>%
    setNames(colnames(snp_counts[, -1])) %>%
    as_tibble %>%
    mutate_all(~ .x / ifelse(!is.na(normalize), .[str_which(snp_counts$snp, normalize)], 1)) %>% # added to normalize by C-T proportions
    add_column(snp = snp_counts$snp, .before = 1)

  snp_props
}


#' Set dimensions of a plot in a Jupyter notebook.
set_dim <- function(width, height, res = 300)
{
  options(
    repr.plot.width = width,
    repr.plot.height = height,
    repr.plot.res = res
  )
}


# df <- tibble(REF = c("A", "G", "T"), pileup = c("cccgt", ",.,,,a", "..,$gc"))
# df %>% mutate(p = count_bases(pileup, REF)) %>% unnest

#' Count bases in a vector of pileup strings.
#' @import purrr
count_bases <- function(pileup, ref) {
  pileup <- tolower(pileup) %>% strsplit("")
  map2(pileup, ref, function(x, y) {
    counts <- factor(x, levels = c("a", "c", "g", "t", ".", ",")) %>%
      table %>%
      as.data.frame %>%
      setNames(c("base", "count")) %>%
      spread(base, count)
    counts$ref_counts <- 0
    if ("." %in% colnames(counts)) {
      counts$ref_counts <- counts[["."]]
      counts[["."]] <- NULL
    }
    if ("," %in% colnames(counts)) {
      counts$ref_counts <- counts$ref_counts + counts[[","]]
      counts[[","]] <- NULL
    }
    counts[[tolower(y)]] <- counts$ref_counts
    counts %>%
      select(-ref_counts)
  })
}
