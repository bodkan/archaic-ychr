#!/usr/bin/env bash

# This scripts takes a subset of sites in given B-team Y chromosome VCF
# that fall in target capture regions. Importantly, since the Y chromosome
# as been called as diploid, it takes only sites which are homozygous, making
# them haploid in the process.

bteam_path=/mnt/454/HighCovNeandertalGenome/1_Extended_VCF

sample_id=$1
info_table=$2
targets_bed=$3

sample_pop=`grep $sample_id $info_table | cut -f3`

bcftools view -g hom -V indels ${bteam_path}/${sample_id}/${sample_id}.hg19_1000g.Y.mod.vcf.gz -R $targets_bed \
    | grep -v "LowQual" \
    | bcftools reheader -s <(echo -e $sample_pop | cat) \
    | sed 's/0\/0/0/; s/1\/1/1/' \
    | bcftools annotate -x INFO,FORMAT -Oz -o vcf/bteam_${sample_id}.vcf.gz
