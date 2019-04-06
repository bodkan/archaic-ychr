#' Assign mutations to all possible branches of a tree with 4 leaves.
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


#' Count the occurences of different site patterns.
sum_patterns <- function(gt, w, x, y, z) {
    df <- gt[, c(w, x, y, z)] %>% .[complete.cases(.), ]
    n <- if (gt$chrom[1] == "simY") max(gt$pos) else nrow(df)
    site_patterns(df, w, x, y, z) %>% summarise_all(sum) %>% mutate(total = n)
}


#' Calculate T_MRCA of all pairs of Africans and non-Africans, using a set
#' of EMHs of known age for branch shortening-based mutation rate estimation.
calculate_tafr <- function(gt, samples) {
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
        tmrca_afr = (tmrca_ad + tmrca_f) / 2
    ) %>%
    select(afr, ref, tmrca_afr, tmrca_ad, tmrca_f, mut_rate, everything())
}


#' Calculate T_MRCA of each of the specified archaic individuals
#' with modern human Y chromosomes.
calculate_tarch <- function(gt, samples, tafr) {
    refs <- filter(samples, pop != "Africa", pop != "EMH")$name
    afrs <- filter(samples, pop == "Africa")$name
    archaics <- colnames(select(gt, -c(chrom, pos, REF, ALT, chimp), -one_of(samples$name)))

    tafr_ui <- filter(tafr, emh == "ustishim") %>%
        select(afr, ref, mut_rate, tmrca_ad, tmrca_f, tmrca_afr)

    map_dfr(refs, function(ref) {
        map_dfr(afrs, function(afr) {
            map_dfr(archaics, function(arch) {
                sum_patterns(gt, w = "chimp", x = arch, y = ref, z = afr) %>%
                    mutate(arch = arch, afr = afr, ref = ref)
            })
        })
    }) %>%
    inner_join(tafr_ui, by = c("afr", "ref")) %>%
    mutate(
      p = a / (a + d + e),
      alpha = (1 + p) / (1 - p),
      tmrca_arch = tmrca_afr * alpha
    ) %>%
    select(arch, afr, ref, tmrca_arch, alpha, mut_rate, tmrca_afr, tmrca_ad, tmrca_f, everything())
}


#' Calculate African and archaic TMRCAs.
calculate_tmrca <- function(gt) {
    if ("simY" %in% gt$chrom)
      samples <- read_siminfo(gt)
    else
      samples <- read_info(gt)

    tafr <- calculate_tafr(gt, samples)
    tarch <- calculate_tarch(gt, samples, tafr)

    tarch
}

