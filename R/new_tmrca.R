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

  map_dfr(refs, ~ sum_patterns(gt, w = "chimp", x = "a00", y = .x, z = "ustishim") %>%
            mutate(ref = .x)) %>%
    select(ref, everything()) %>%
    add_mut_rate %>%
    add_tafr
}

run_step2 <- function(gt, step1) {
  tafr <- select(step1, ref, mut_rate, tmrca_afr)

  samples <- if ("simY" %in% gt$chrom) read_siminfo(gt) else read_info(gt)

  refs <- filter(samples, pop != "Africa", pop != "EMH")$name
  archaic <- colnames(select(gt, -c(chrom, pos, REF, ALT, chimp), -one_of(samples$name)))

  map_dfr(refs, ~ sum_patterns(gt, w = "chimp", x = archaic, y = .x, z = "a00") %>%
            mutate(arch = archaic, ref = .x)) %>%
    inner_join(tafr, by = "ref") %>%
    add_tarch_mendez %>%
    add_tarch_new %>%
    select(arch, ref, everything())
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


add_tafr_ci <- function(df) {
  simulate_counts(df) %>%
    add_mut_rate %>%
    add_tafr %>%
    select(tmrca_afr, mut_rate) %>%
    calculate_ci %>%
    bind_cols(df, .) %>%
    select(starts_with("tmrca"), starts_with("mut_rate"), everything())
}


add_tarch_ci <- function(df, tafr) {
  simulate_counts(df) %>%
    mutate(tmrca_afr = tafr) %>%
    add_tarch_mendez %>%
    add_tarch_new %>%
    select(tmrca_mendez, tmrca_new) %>%
    calculate_ci %>%
    bind_cols(df, .) %>%
    select(starts_with("tmrca"), starts_with("mut_rate"), everything())
}


calculate_ci <- function(df) {
  df %>%
    gather(stat, value) %>%
    group_by(stat) %>%
    summarise(
      low = quantile(value, 0.025, na.rm = TRUE),
      high = quantile(value, 0.975, na.rm = TRUE)
    ) %>%
    gather(boundary, value, -stat) %>%
    unite(temp, stat, boundary) %>%
    spread(temp, value)
}
