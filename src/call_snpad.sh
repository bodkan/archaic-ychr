#!/usr/bin/env bash

# Modified from Janet's snpAD README to be run not just on Mez2,
# but also on merged Den4+Den8 data.

name=$1
bam=$2
vcf=$3

if [[ $# != 3 ]]; then
    echo "Insufficient number of arguments"
    exit
fi

prefix=`basename $bam`

dir=tmp/snpad/${name}
mkdir -p $dir

# prepare snpAD input files from bam 
/home/cesare_filippo/bin/Bam2snpAD \
    -r Y \
    -f /mnt/solexa/Genomes/hg19_evan/whole_genome.fa \
    -Q 25 $bam \
    > ${dir}/${prefix}.input.snpad

/home/cesare_filippo/src/intersectbed.pl \
    ${dir}/${prefix}.input.snpad \
    <(tabix /mnt/454/HCNDCAM/Hengs_Alignability_Filter/hs37m_filt35_50.bed.gz Y) \
    > ${dir}/${prefix}.input.tab

# generate priors and errors
/home/pruefer/bin/snpAD \
    -c 25 \
    -o ${dir}/${prefix}.priors.txt \
    -O ${dir}/${prefix}.errors.txt \
    ${dir}/${prefix}.input.tab \
    > ${dir}/${prefix}.prior_errors.log 2>&1

# call
/home/pruefer/bin/snpADCall \
    -N $name \
    -e ${dir}/${prefix}.errors.txt \
    -p "`cat ${dir}/${prefix}.priors.txt`" \
    ${dir}/$prefix.input.tab \
    | bgzip \
    > $vcf
