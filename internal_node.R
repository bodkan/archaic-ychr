library(tidyverse)
library(magrittr)
library(here)
library(furrr)
library(binom)

devtools::load_all(".")

plan(multiprocess)

options(future.globals.maxSize = 1500*1024^2)


highcov_all %<-% read_genotypes("mez2", "full", mindp = 3, maxdp = 0.98, nodmg = F, var_only = T)
highcov_nodmg %<-% read_genotypes("mez2", "full", mindp = 3, maxdp = 0.98, nodmg = T, var_only = T)


nrow(highcov_all)
nrow(highcov_nodmg)


date_split <- function(ind_a, ind_b, gt, ancestral_ind, mindp = 3) {
    ind_a <- str_replace(ind_a, "-", "_")

    ## the individual of interest is not in the highcov VCF => load archaic VCF
    if (!ind_a %in% colnames(gt)) {
      capture <- ifelse(ind_a == "elsidron2", "lippold", "full")
      ind_gt <- read_vcf(here(paste0("data/vcf/", capture, "_", ind_a, ".vcf.gz")), mindp = mindp, maxdp = 0.98) %>% filter(!is.na(!!sym(ind_a)))
    } else { # the individual of interest is in the modern human VCF
        ind_gt <- select(gt, chrom, pos, REF, ALT, !!ind_a)
        gt <- select(gt, -!!ind_b)
    }
    
    gt <- gt[, c("chrom", "pos", "REF", "ALT", ind_b, ancestral_ind)]
    ancestral_freq <- gt[, ancestral_ind] %>% rowMeans
    from_derived <- filter(gt, !!sym(ind_b) == 1, ancestral_freq == 0)
    
    joined_gt <- inner_join(from_derived, ind_gt, by = c("chrom", "pos", "REF")) %>%
        filter(complete.cases(.)) %>%
        select(-ALT.y) %>%
        rename(ALT = ALT.x)

    joined_gt
}















date_split2 <- function(ind_a, ind_b, gt, ancestral_ind, mindp = 3) {
    capture <- ifelse(ind_a == "elsidron2", "lippold", "full")

    ind_gt <- here(paste0("data/vcf/", capture, "_", ind_a, ".vcf.gz")) %>%
        read_vcf(mindp = mindp, maxdp = 0.98) %>%
        filter(!is.na(!!sym(ind_a)))

    gt <- gt[, c("chrom", "pos", "REF", "ALT", ind_b, ancestral_ind)]
    ancestral_freq <- gt[, ancestral_ind] %>% rowMeans
    from_derived <- filter(gt, !!sym(ind_b) == 1, ancestral_freq == 0)

    joined_gt <- inner_join(from_derived, ind_gt, by = c("chrom", "pos", "REF")) %>%
        filter(complete.cases(.)) %>%
        select(-ALT.y) %>%
        rename(ALT = ALT.x)

    joined_gt
}




date_split("spy1", "mez2", highcov_all, c("chimp", "a00", "S_French_1"), mindp = 3) %>%
    group_by(spy1) %>% tally %>% mutate(prop = n / sum(n))
date_split2("spy1", "mez2", highcov_nodmg, c("chimp", "a00", "S_French_1"), mindp = 3) %>%
    group_by(spy1) %>% tally %>% mutate(prop = n / sum(n))


date_split2("spy1", "mez2", highcov_all, c("chimp", "a00", "S_French_1"), mindp = 2) %>%
    group_by(spy1) %>% tally %>% mutate(prop = n / sum(n))
date_split2("spy1", "mez2", highcov_nodmg, c("chimp", "a00", "S_French_1"), mindp = 2) %>%
    group_by(spy1) %>% tally %>% mutate(prop = n / sum(n))


date_split2("spy1", "mez2", highcov_all, c("chimp", "a00", "S_French_1"), mindp = 1) %>%
    group_by(spy1) %>% tally %>% mutate(prop = n / sum(n))
date_split2("spy1", "mez2", highcov_nodmg, c("chimp", "a00", "S_French_1"), mindp = 1) %>%
    group_by(spy1) %>% tally %>% mutate(prop = n / sum(n))


date_split2("elsidron2", "mez2", highcov_all, c("chimp", "a00", "S_French_1"), mindp = 3) %>%
    group_by(elsidron2) %>% tally %>% mutate(prop = n / sum(n))
date_split2("elsidron2", "mez2", highcov_nodmg, c("chimp", "a00", "S_French_1"), mindp = 3) %>%
    group_by(elsidron2) %>% tally %>% mutate(prop = n / sum(n))





