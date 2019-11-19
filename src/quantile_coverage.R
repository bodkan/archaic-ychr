args <- commandArgs(trailingOnly = T)
if (length(args) != 2)
    stop("Path to a VCF file and a quantile value required as arguments", call. = F)

path <- args[1]
cutoff <- as.numeric(args[2])

if (!file.exists(path))
    stop(paste0("VCF file '", path, "' does not exist"))

if (cutoff <= 0 | cutoff >= 1)
    stop("Quantile cutoff must lie in the interval (0, 1)")

suppressPackageStartupMessages({
library(VariantAnnotation)
})

vcf <- readVcf(path, row.names = FALSE)

dp <- geno(vcf)$DP

quant <- apply(dp, 2, function(i) quantile(i, cutoff, na.rm = TRUE))

cat(quant)













