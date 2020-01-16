# Read and filter a single VCF file.
read_vcf <- function(path, mindp, maxdp = 0.98, var_only = FALSE, nodmg = FALSE, bed_filter = NA) {
  vcf <- VariantAnnotation::readVcf(path, row.names = FALSE)
  gr <- GenomicRanges::granges(vcf)

  if (!"DP" %in% names(VariantAnnotation::geno(vcf))) {
    dims <- dim(VariantAnnotation::geno(vcf)$GT)
    mask <- matrix(TRUE, dims[1], dims[2])
  } else {
    dp <- VariantAnnotation::geno(vcf)$DP
    # apply min and max coverage filters
    mask <- apply(dp, 2, function(i) ifelse(i >= mindp & i <= quantile(i, maxdp, na.rm = TRUE), TRUE, FALSE))
    if ("chimp" %in% colnames(mask)) mask[, "chimp"] <- TRUE
  }

  biallelic_pos <- IRanges::elementNROWS(gr$ALT) == 1

  gt <- VariantAnnotation::geno(vcf)$GT
  # process snpAD diploid calls correctly
  gt <- gt %>% replace(. == "0/0", "0") %>% replace(. == "1/1", "1") %>% replace(nchar(.) > 1, ".")
  # keep genotypes only for sites that are present, or pass the filtering
  gt <- gt %>% replace(. == ".", NA) %>% replace(!mask, NA)
  mode(gt) <- "numeric"

  gt_df <- tibble::as_tibble(gt) %>% filter(biallelic_pos)
  info_df <- tibble::tibble(
    chrom = as.character(GenomicRanges::seqnames(gr))[biallelic_pos],
    pos = GenomicRanges::start(gr)[biallelic_pos],
    REF = as.character(gr$REF)[biallelic_pos],
    ALT = as.character(unlist(gr$ALT[biallelic_pos, ]))
  )

  df <- dplyr::bind_cols(info_df, gt_df) %>% dplyr::filter(REF != "N")

  # sanitize the sample names
  colnames(df) <- str_replace_all(colnames(df), "-", "_")

  # not_missing <- select(df, -chrom:-ALT) %>% is.na %>% apply(MARGIN = 1, all)

  # df <- df %>% filter(!not_missing)

  if (var_only) df <- filter(df, ALT != "")

  if (nodmg) {
    df <- filter(df, ALT == "" | !((REF == "C" & ALT == "T") | (REF == "G" & ALT == "A")))
  }

  if (!is.na(bed_filter)) {
      bed_filter <- rtracklayer::import.bed(bed_filter)
    df <- df %>%
      GenomicRanges::makeGRangesFromDataFrame(
        start.field = "pos",
        end.field = "pos",
        keep.extra.columns = TRUE
      ) %>%
      IRanges::subsetByOverlaps(bed_filter) %>%
      as.data.frame %>%
      dplyr::select(-end, -width, -strand) %>%
      dplyr::rename(chrom = seqnames, pos = start) %>%
      tibble::as_tibble()
  }

  df
}


#' Read genotypes from a VCF file, returning a data frame object.
#' @param capture Capture set (full, lippold, exome).
#' @param archaic Path to a low-coverage archaic Y chromosome VCF.
#' @param mindp Minimum coverage at each site.
#' @param maxdp Maximum coverage at each site (specified as a proportion of an
#'   upper tail of the entire coverage distribution).
#' @import stringr dplyr purrr tibble
read_genotypes <- function(archaic, capture, mindp, maxdp = 0.98, var_only = FALSE, nodmg = FALSE, bed_filter = NA) {
  archaic_vcf <- here::here(paste0("data/vcf/", capture, "_", archaic, ".vcf.gz"))
  highcov_vcf <- here::here(paste0("data/vcf/", ifelse(capture == "test", "full", capture), "_modern.vcf.gz"))

  archaic_df <- read_vcf(archaic_vcf, mindp = mindp, maxdp = maxdp, bed_filter = bed_filter)
  highcov_df <- read_vcf(highcov_vcf, mindp = mindp, maxdp = maxdp, bed_filter = bed_filter)

  df <- dplyr::right_join(archaic_df, highcov_df, by = c("chrom", "pos" ,"REF"), suffix = c("_arch", "_modern"))

  # remove third alleles from the archaic human sample
  archaic_name <- archaic_df %>% colnames %>% .[length(.)]
  df <- mutate(df, !!archaic_name := ifelse((ALT_modern != "" & ALT_arch != "" & ALT_modern != ALT_arch), NA, df[[archaic_name]]))

  # collapse ALT columns discovered in modern and archaic samples
  df <- mutate(df, ALT = case_when(
      ALT_modern != ""                    ~ ALT_modern,
      ALT_modern == "" & !is.na(ALT_arch) ~ ALT_arch,
      TRUE                                ~ "")
    ) %>%
    select(chrom, pos, REF, ALT, everything()) %>%
    select(-ALT_modern, -ALT_arch)

  if (var_only) df <- filter(df, ALT != "")

  if (nodmg) {
    df <- filter(df, ALT == "" | !((REF == "C" & ALT == "T") | (REF == "G" & ALT == "A")))
  }

  df
}


#' Read table of genotypes simulated by msprime.
#' @import dplyr readr
read_simgt <- function(path) {
  suppressMessages(read_tsv(path)) %>%
    mutate(chrom = "simY", pos = round(pos), REF = "A", ALT = "T") %>%
    select(chrom, pos, REF, ALT, everything())
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
