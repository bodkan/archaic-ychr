#!/bin/bash

mkdir revisions
cd revisions

conda activate archaic-ychr

# restrict to deaminate reads only
for bam in den4 den8 den_merged spy1 mez2; do
    ~mmeyer/perlscripts/solexa/analysis/filterBAM.pl -p5 0,1,2 -p3 0,-1,-2 -suffix deam ../data/bam/full_${bam}.bam
    mv full_${bam}.deam.bam full_${bam}_deam.bam
done

# determine patterns of deamination as a sanity check
for bam in den4 den8 den_merged spy1 mez2; do
    ~mmeyer/perlscripts/solexa/analysis/substitution_patterns.pl full_${bam}_deam.bam
done

# index
for bam in den4 den8 den_merged spy1 mez2; do
    samtools index full_${bam}_deam.bam
done

# run consensus genotype calling
for bam in den4 den8 den_merged spy1 mez2; do
    python3 /mnt/expressions/mp/bam-caller/bam-caller.py --bam full_${bam}_deam.bam \
         --strategy majority --proportion 0.9 --mincov 1 --minbq 20 --minmq 25 \
         --sample-name ${bam}_deam --output full_${bam}_deam
    bgzip full_${bam}_deam.vcf
    tabix full_${bam}_deam.vcf.gz
done

# copy the new genotype calls to the main VCF directory
cp *.vcf.gz* ../data/vcf

# download the Davalos regions
curl https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-018-4945-x/MediaObjects/12864_2018_4945_MOESM1_ESM.bed \
    | sed 's/chr//' \
    > davalos.bed

bedtools intersect -a davalos.bed -b /mnt/454/HCNDCAM/Hengs_Alignability_Filter/hs37m_filt35_50.bed.gz \
    > davalos_map50.bed
