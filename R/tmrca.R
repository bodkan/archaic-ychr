add_mut_rate <- function(counts) {
  mutate(counts, mut_rate = (d - e) / 45000 / total)
}

add_tafr <- function(counts) {
  mutate(counts, tmrca_afr = f / (total * mut_rate))
}

add_tarch_mendez <- function(counts) {
  mutate(counts, tmrca_mendez = tmrca_afr * (2 * a + d + e) / (d + e))
}

add_tarch_new <- function(counts) {
  mutate(counts, tmrca_new = tmrca_afr * (a + d) / d)
}


run_step1 <- function(gt) {
  samples <- if ("simY" %in% gt$chrom) read_siminfo(gt) else read_info(gt)

  refs <- filter(samples, pop != "Africa", pop != "EMH")$name
  afrs <- filter(samples, pop == "Africa")$name

  map_dfr(afrs, function(afr) {
  map_dfr(refs, ~ sum_patterns(gt, w = "chimp", x = afr, y = .x, z = "ustishim") %>%
            mutate(afr = afr, ref = .x))
  }) %>%
    select(afr, ref, everything()) %>%
    add_mut_rate %>%
    add_tafr
}

run_step2 <- function(gt, step1) {
  tafr <- select(step1, afr, ref, mut_rate, tmrca_afr)

  samples <- if ("simY" %in% gt$chrom) read_siminfo(gt) else read_info(gt)

  refs <- filter(samples, pop != "Africa", pop != "EMH")$name
  archaic <- colnames(select(gt,
                             -c(chrom, pos, REF, ALT, chimp),
                             -one_of(samples$name),
                             -starts_with("S_")))
  afrs <- filter(samples, pop == "Africa")$name

  map_dfr(afrs, function(afr) {
  map_dfr(refs, ~ sum_patterns(gt, w = "chimp", x = archaic, y = .x, z = afr) %>%
            mutate(arch = archaic, afr = afr, ref = .x))
  }) %>%
    inner_join(tafr, by = c("afr", "ref")) %>%
    add_tarch_mendez %>%
    add_tarch_new %>%
    select(arch, afr, ref, everything())
}


simulate_counts <- function(counts, n = 1000) {
  tibble(
    a = rpois(n, counts$a),
    b = rpois(n, counts$b),
    c = rpois(n, counts$c),
    d = rpois(n, counts$d),
    e = rpois(n, counts$e),
    f = rpois(n, counts$f),
    total = counts$total
  )
}


#' Take A00 TMRCA and mutation rate estimates and simulate
#' confidence intervals by bootstraping for both of them.
#' @import tidyr
add_step1_ci <- function(df, per_ind) {
  if (!per_ind) {
    columns <- "afr"
    df <- group_by(df, afr) %>% summarise_if(is.numeric, mean)
  } else {
    columns <- c("afr", "ref")
  }

  ci <- df %>%
    nest(counts = c(a, b, c, d, e, f, total)) %>%
    mutate(sims = map(counts, simulate_counts)) %>%
    unnest(sims) %>%
    add_mut_rate %>%
    add_tafr %>%
    select(one_of(columns), tmrca_afr, mut_rate) %>%
    calculate_ci(per_ind)

  inner_join(df, ci, by = columns) %>%
    select(one_of(columns), starts_with("tmrca"), starts_with("mut_rate"))
}

#' Take archaic Y chromosome TMRCA estimates and simulate
#' confidence intervals by bootstraping of branch counts.
add_step2_ci <- function(df, per_ind) {
  columns <- if (per_ind) c("arch", "afr", "ref") else c("arch", "afr")

  ci <- df %>%
    nest(counts = c(a, b, c, d, e, f, total)) %>%
    mutate(sims = map(counts, simulate_counts)) %>%
    unnest(sims) %>%
    add_tarch_mendez %>%
    add_tarch_new %>%
    select(one_of(columns), tmrca_mendez, tmrca_new) %>%
    calculate_ci(per_ind)

  if (!per_ind) {
    df <- select(df, afr, arch, matches("mendez|new")) %>%
      group_by_at(vars(one_of(columns))) %>%
      summarise_if(is.numeric, mean)
  }

  inner_join(df, ci, by = columns) %>%
    select(one_of(columns), starts_with("tmrca"))
}


#' Calculate empirical confidence intervals from bootstrapped
#' values of simulated branch counts.
calculate_ci <- function(df, per_ind) {
  columns <- c("afr")
  if ("arch" %in% names(df)) columns <- c(columns, "arch")
  if (per_ind) {
    columns <- c(columns, "ref")
  }

  df %>%
    gather(stat, value, -one_of(columns)) %>%
    group_by_at(vars(one_of(c("stat", columns)))) %>%
    summarise(
      low = quantile(value, 0.025, na.rm = TRUE),
      high = quantile(value, 0.975, na.rm = TRUE)
    ) %>%
    gather(boundary, value, -one_of(c("stat", columns))) %>%
    unite(temp, stat, boundary) %>%
    spread(temp, value)
}

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

