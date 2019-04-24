# Extract VCF info columns from a GRanges object.
vcf_info <- function(gr) {
  tibble::tibble(
    chrom = as.character(GenomicRanges::seqnames(gr)),
    pos = GenomicRanges::start(gr),
    REF = as.character(gr$REF),
    ALT = as.character(unlist(gr$ALT))
  )
}


# arch_path = "data/vcf/lippold_den8.vcf.gz"
# modern_path = "data/vcf/merged_lippold.vcf.gz"
# mindp = 1
# maxdp = 0.975
# var_only = F

#' Read genotypes from a VCF file, returning a data frame object.
#' @param mindp Minimum coverage at each site.
#' @param maxdp Maximum coverage at each site (specified as a proportion of an
#'   upper tail of the entire coverage distribution).
#' @import stringr dplyr purrr tibble
read_gt <- function(arch_path, modern_path, mindp, maxdp = 0.975, var_only = FALSE, tv_only = FALSE) {
  vcf_arch <- VariantAnnotation::readVcf(arch_path)
  vcf_modern <- VariantAnnotation::readVcf(modern_path)

  gr_arch <- GenomicRanges::granges(vcf_arch)
  gr_modern <- GenomicRanges::granges(vcf_modern)

  # read DP information for all samples in both VCFs
  dp_arch <- VariantAnnotation::geno(vcf_arch)$DP
  dp_modern <- VariantAnnotation::geno(vcf_modern)$DP

  # apply min and max coverage filters
  mask_arch <- apply(dp_arch, 2, function(i) ifelse(i >= mindp & i <= quantile(i, maxdp, na.rm = TRUE), TRUE, FALSE))
  mask_modern <- apply(dp_modern, 2, function(i) ifelse(i >= 4 & i <= quantile(i, maxdp, na.rm = TRUE), TRUE, FALSE))
  if ("chimp" %in% colnames(mask_modern)) mask_modern[, "chimp"] <- TRUE

  # keep genotypes only for sites that are present, or pass the filtering
  gt_modern <- VariantAnnotation::geno(vcf_modern)$GT %>% replace(. == ".", NA) %>% replace(!mask_modern, NA)
  gt_arch <- VariantAnnotation::geno(vcf_arch)$GT %>% replace(. == ".", NA) %>% replace(!mask_arch, NA)

  mode(gt_modern) <- "numeric"
  mode(gt_arch) <- "numeric"

  gt_modern <- tibble::as_tibble(gt_modern) %>% dplyr::mutate(reference = 0)
  gt_arch <- tibble::as_tibble(gt_arch)

  # REMOVE SELECT
  df_modern <- vcf_info(gr_modern) %>% dplyr::bind_cols(gt_modern)
  df_arch <- vcf_info(gr_arch) %>% dplyr::bind_cols(gt_arch)

  df <- dplyr::right_join(df_arch, df_modern, by = c("chrom", "pos" ,"REF"), suffix = c("_modern", "_arch"))
  # sanitize the sample names
  colnames(df) <- str_replace_all(colnames(df), "-", "_")

  # remove tri-allelic sites
  df <- filter(df, !(ALT_modern != "" & ALT_arch != "" & ALT_modern != ALT_arch))

  # collapse ALT columns discovered in modern and archaic samples
  df <- mutate(df, ALT = ifelse(ALT_modern != "", ALT_modern, ALT_arch)) %>%
    select(chrom, pos, REF, ALT, everything()) %>%
    select(-ALT_modern, -ALT_arch)

  if (var_only) df <- filter(df, ALT != "")

  if (tv_only) {
    df <- filter(res_df,  !((REF == "C" & ALT == "T") | (REF == "G" & ALT == "A")))
  }

  df
}


#' Read table of genotypes simulated by msprime.
#' @import dplyr readr
read_simgt <- function(path) {
  suppressMessages(read_tsv(path)) %>%
    mutate(chrom = "simY", pos = round(pos), REF = "A", ALT = "T") %>%
    select(chrom, pos, REF, ALT, chimp = chimp0, everything())
}


#' Add a given proportion of sequencing/damage errors to a set of samples.
#' samples <- c("arch0", "arch1")
#' @import purrr
add_errors <- function(gt, rate, samples) {
  gt[, samples] <- map_dfr(gt[, samples], function(alleles) {
    mut_pos <- sort(sample(seq_along(alleles), size = length(alleles) * rate))
    alleles[mut_pos] <- as.integer(!as.logical(alleles[mut_pos]))
    alleles
  })
  gt
}
