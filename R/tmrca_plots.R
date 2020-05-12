#' Filter sites/depth/capture regions/individuals
filter_tmrca <- function(stat, afr, mindp, sites, filt, exclude = NA) {
    filter(
        tmrca_df,
        afr == !!afr,
        filt == !!filt,
        sites == !!sites & dp == mindp & capture == "full" & (str_detect(arch, "mez2_dp") | arch %in% c("den4", "den8", "mez2", "spy1", "shotgun_spy1", "shotgun_mez2", "mez2_snpad", "den_snpad", "den", "spy1_snpad", "den4_snpad", "den8_snpad") | str_detect(arch, "deam")) |
        sites == !!sites & dp == mindp & capture == "lippold" & arch == "elsidron2" |
        dp %in% c(1, 3) & sites == "all" & capture == "exome" & arch == "elsidron1",
        !arch %in% exclude
    ) %>%
        mutate(arch = case_when(arch == "elsidron1" & dp == 1 ~ "elsidron_dp1",
                                arch == "elsidron1" & dp == 3 ~ "elsidron_dp3",
                                TRUE ~ arch))
}



#' Combine archaic and A00 TMRCA data frames for further analysis
#' Calculate confidence intervals for each reference individual by bootstrapping
get_ind_tmrca <- function(stat, afr, mindp, sites, filt, exclude = NA) {
    arch_ind_ci <- filter_tmrca(stat, afr, mindp, sites, filt, exclude) %>%
        add_step2_ci(per_ind = TRUE) %>%
        select(name = arch, ref, contains(stat)) %>%
        setNames(c("name", "ref", "tmrca", "high", "low"))

    afr_ind_ci <- step1 %>%
        filter(afr == !!afr) %>%
        add_step1_ci(per_ind = TRUE) %>%
        select(name = afr, ref, contains("tmrca_afr")) %>%
        filter(!name %in% exclude) %>%
        setNames(c("name", "ref", "tmrca", "high", "low"))

    bind_rows(arch_ind_ci, afr_ind_ci) %>%
        assign_set %>%
        mutate(set = factor(set, levels = c("Denisovan", "Neanderthal", "African")))
}



get_overall_tmrca <- function(stat, afr, mindp, sites, filt, exclude = NA) {
    arch_overall <- filter_tmrca(stat, afr, mindp, sites, filt, exclude) %>%
        add_step2_ci(per_ind = FALSE) %>%
        select(name = arch, contains(stat)) %>%
        setNames(c("name", "tmrca", "high", "low"))
   
    # overall A00 TMRCA values and CIs
    afr_overall <- step1 %>%
        filter(afr == !!afr) %>%
        add_step1_ci(per_ind = FALSE) %>%
        select(name = afr, contains("tmrca_afr")) %>%
        filter(!name %in% exclude) %>%
        select(name, tmrca = tmrca_afr, low = tmrca_afr_low, high = tmrca_afr_high)

    bind_rows(arch_overall, afr_overall) %>%
        assign_set %>%
        mutate(set = factor(set, levels = c("Denisovan", "Neanderthal", "African")))
}



plot_tmrca <- function(stat, afr, mindp, sites, filt, ylabel = TRUE, exclude = NA, title = NA) {
    colors <- scales::hue_pal()(3)[c(2, 3, 1)]

    ind_tmrca <- get_ind_tmrca(stat, afr, mindp, sites, filt, exclude) %>%
        mutate(name = fix_name(name)) %>%
        mutate(name = str_replace(name, " \\(560 kb\\)", "")) %>%
        mutate(name = fct_reorder(name, as.numeric(set))) %>%
        mutate(name = fct_relevel(name,
                                  "El Sidrón 1253 (118 kb, filtered)",
                                  "El Sidrón 1253 (118 kb, unfiltered)", after = Inf)) %>%
        mutate(name = fct_relevel(name, "Denisova 4 & 8 (snpAD)", after = 2))

    overall_tmrca <- get_overall_tmrca(stat, afr, mindp, sites, filt, exclude) %>%
        mutate(name = fix_name(name)) %>%
        mutate(name = str_replace(name, "\\(560 kb\\)", "")) %>% 
        mutate(name = fct_reorder(name, as.numeric(set))) %>%
        mutate(name = fct_relevel(name,
                                  "El Sidrón 1253 (118 kb, filtered)",
                                  "El Sidrón 1253 (118 kb, unfiltered)", after = Inf)) %>%
        mutate(name = fct_relevel(name, "Denisova 4 & 8 (snpAD)", after = 2))

    p <- ggplot() +
        # individual CIs
        geom_linerange(data = ind_tmrca,
                       aes(name, ymin = low, ymax = high, group = ref, color = set),
                       position = position_dodge(width = 0.7), size = 1/5, alpha = 1.0) +

        # individual points
        geom_point(data = ind_tmrca, aes(name, tmrca, group = ref, color = set),
                      position = position_dodge(width = 0.7), size = 1/5, alpha = 1.0) +
        # average
        geom_segment(data = overall_tmrca, size = 1/3, linetype = 2,
                     aes(x = as.numeric(name) - 0.42, y = tmrca, xend = as.numeric(name) + 0.42, yend = tmrca)) +
        # upper bar
        geom_segment(data = overall_tmrca, size = 1/4,
                     aes(x = as.numeric(name) - 0.35, y = high, xend = as.numeric(name) + 0.35, yend = high)) +
        # lower bar
        geom_segment(data = overall_tmrca, size = 1/4,
                     aes(x = as.numeric(name) - 0.35, y = low, xend = as.numeric(name) + 0.35, yend = low)) +
        theme(axis.text.x = element_text(hjust = 1, angle = 45)) + 
        theme_classic() +
        scale_y_continuous(labels = comma) +
        theme(
            legend.position = "none",
            axis.text.x = element_text(hjust = 1, angle = 30, size = 13),
            axis.title.x = element_blank(),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 13),
            plot.subtitle = element_text(hjust = 1)
        ) +
        background_grid(major = "xy", minor = "xy", size.major = 0.2, size.minor = 0.2,
                        color.major = rgb(0.8, 0.8, 0.8, alpha = 0.5),
                        color.minor = rgb(0.8, 0.8, 0.8, alpha = 0.5)) +
        coord_capped_cart(ylim = c(0, 1e6), left = "both") +
    
        scale_color_manual(values = colors)

    if (ylabel) {
        p <- p + ylab("TMRCA with non-Africans [years ago]")
    } else {
        p <- p + theme(axis.title.y = element_blank())
    }

    if (is.na(title)) return(p)

    p + ggtitle("", subtitle = title)
}
