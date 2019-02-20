#!/usr/bin/env bash

if [[ $# -ne 2 ]]; then
    echo -e "Usage:\n\t./chimp_vcf.sh <positions file> <output vcf>"
    exit 1;
fi

positions=$1
output=$2

echo "##fileformat=VCFv4.1" >> $output
echo "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >> $output
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tChimp" >> $output

for chrom in {1..22} X Y; do
    grep -w "^${chrom}" $positions \
        | /home/pruefer/bin/BamSNPAddMaf \
                /mnt/sequencedb/ucsc/goldenPath/hg19/vsPanTro4/axtNet/mafnet/chr${chrom}.maf \
                hg19 panTro4 \
        | awk '{print toupper($0)}' \
        | awk -vOFS="\t" '
            {
                if ($3 == "N") {
                    next
                }
                else if ($4 == "N" || $4 == "-") {
                    alt = "."
                    gt = "."
                } else if ($3 == $4) {
                    alt = "."
                    gt = "0"
                } else {
                    alt = $4
                    gt = "1"
                }
                { print $1, $2, ".", $3, alt, ".", ".", ".", "GT", gt}
            }
        '
done >> $output
