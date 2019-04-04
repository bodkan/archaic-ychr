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
