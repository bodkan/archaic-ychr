#!/usr/bin/env bash

if [[ $# -ne 1 ]]; then
    echo "Missing path to a VCF file"
    exit
fi

vcf=$1

cutoff=`Rscript src/quantile_coverage.R ${vcf} 0.98`

echo "Filtering ${vcf} using DP >= 3 && DP <= ${cutoff}..."

bcftools view -i "DP >= 3 && DP <= ${cutoff}" $vcf -Oz -o tmp/vcf_fasta/`basename $vcf`; tabix tmp/vcf_fasta/`basename $vcf`
