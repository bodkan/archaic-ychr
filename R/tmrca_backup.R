add_tafr_ci <- function(df) {
  simulate_counts(df) %>%
    add_mut_rate %>%
    add_tafr %>%
    select(tmrca_afr, mut_rate) %>%
    calculate_ci2 %>%
    bind_cols(df, .) %>%
    select(starts_with("tmrca"), starts_with("mut_rate"), everything())
}



add_tarch_ci <- function(df, tafr) {
  simulate_counts(df) %>%
    mutate(tmrca_afr = tafr) %>%
    add_tarch_mendez %>%
    add_tarch_new %>%
    select(tmrca_mendez, tmrca_new) %>%
    calculate_ci2 %>%
    bind_cols(df, .) %>%
    select(starts_with("tmrca"), starts_with("mut_rate"), everything())
}


calculate_ci2 <- function(df) {
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
