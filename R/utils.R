#' Generate table of sample names and their population assignments and ages
#' from the simulated data.
read_siminfo <- function(gt) {
  emh <- tibble(name = str_subset(colnames(gt), "ustishim"), age = 45000, pop = "EMH")
  afr <- tibble(name = str_subset(colnames(gt), "afr"), age = 0, pop = "Africa")
  westeur <- tibble(name = str_subset(colnames(gt), "eur"), age = 0, pop = "WestEur")
  bind_rows(emh, afr, westeur)
}


#' Generate table of sample names and their population assignments and ages.
read_info <- function(gt) {
  emh <- tibble(
    name = c("ustishim"), #"bichon", "kk1", "loschbour", "mota"),
    age = c(45000), #, 13665, 9720, 8050, 4500),
    pop = "EMH"
  )

  modern <- tibble(
    name = c(str_subset(colnames(gt), "^S_"), "a00"), #, "a00_1", "a00_2"),
    pop = case_when(
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
    mutate(snp = paste0(REF, "-", ALT)) %>%
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


#' Read coverage info for a given BAM file in a defined set of regions.
get_coverage <- function(bed, bam) {
  paste("bedtools coverage -a", bed, "-b", bam, "-d") %>%
    pipe %>%
    read_tsv(col_names = c("chr", "start", "end", "pos", "coverage"),
             col_types = "ciiii") %>%
    makeGRangesFromDataFrame(starts.in.df.are.0based = TRUE,
                             keep.extra.columns = TRUE)
}


#' Change sample name to its full form.
fix_name <- function(name, coverage = FALSE) {
  new_name <- case_when(
    name == "den4" ~ "Denisova 4",
    name == "den8" ~ "Denisova 8",
    name == "spy1" ~ "Spy 94a",
    name == "mez2" ~ "Mezmaiskaya 2",
    name == "mez2_snpad" ~ "Mezmaiskaya 2 (snpAD)",
    name == "shotgun_mez2" ~ "Mezmaiskaya 2 (shotgun)",
    name == "shotgun_spy1" ~ "Spy 94a (shotgun)",
    name == "elsidron1" ~ "El Sidrón 1253 (118 kb)",
    name == "elsidron_dp1" ~ "El Sidrón 1253 (118 kb, unfiltered)",
    name == "elsidron_dp3" ~ "El Sidrón 1253 (118 kb, filtered)",
    name == "elsidron2" ~ "El Sidrón 1253",
    name == "a00" ~ "A00 lineage",
    name == "ustishim" ~ "Ust'-Ishim",
    TRUE ~ name
  )

  cov_df <- tribble(
    ~name, ~coverage,
    "Spy 94a", 0.8253635,
    "Denisova 4", 1.5451500,
    "El Sidrón 1253 (118 kb)",	3.1318824,
    "Denisova 8", 3.4819494,
    "El Sidrón 1253", 7.8449213,
    "Mezmaiskaya 2", 14.2548400
  ) %>% mutate(coverage = round(coverage, 1))

  if (coverage)
    new_name <- map_chr(new_name, ~ paste0(.x, " (", filter(cov_df, name == .x)$coverage, "X)"))

  new_name
}

summean <- function(df) summarise_if(df, is.numeric, mean)

#' Add a column with a sample group.
assign_set <- function(df) {
  ungroup(df) %>%
  mutate(set = case_when(name == "a00" ~ "A00 lineage",
             name %in% c("den4", "den8") ~ "Denisovan",
             name %in% c("spy1", "mez2", "shotgun_mez2", "shotgun_spy1", "mez2_snpad") ~ "Neanderthal",
             str_detect(name, "elsidron") ~ "Neanderthal",
             str_detect(name, "shotgun") ~ "Neanderthal (shotgun)",
             TRUE ~ "other")) %>%
  mutate(name = fix_name(name)) %>%
  mutate(name = fct_relevel(name, "Denisova 4", "Denisova 8", "Spy 94a (shotgun)", "Spy 94a",
                            "El Sidrón 1253", "Mezmaiskaya 2 (shotgun)", "Mezmaiskaya 2",
                            "Mezmaiskaya 2 (snpAD)", "A00"))
}


#' Pipe operator
#'
#' Added via usethis::use_pipe().
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL
