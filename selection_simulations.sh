#!/usr/bin/env bash

for chrom in Y; do
for gene_total in `seq 100000 100000 2000000`; do
for admix_rate in `seq 0.01 0.01 0.1`; do
for admix_time in `seq 200000 50000 450000`; do
for rep in `seq 1 25`; do
    N="modern2neand_sel_chr${chrom}_seq${gene_total}_rate${admix_rate}_time${admix_time}_rep${rep}"
    qsub -V -b yes -cwd -j y -l virtual_free=5G,h_vmem=5G -o tmp/sge/${N}.out -N $N \
        slim \
            -d "chrom='${chrom}'" \
            -d mut_rate=1.85e-08 \
            -d chrom_length=60000000 \
            -d gene_count=100 \
            -d gene_total=${gene_total} \
            -d admix_time=${admix_time} \
            -d modern_neand=${admix_rate} \
            -d neand_modern=0 \
            -d "output='data/newsim/${N}.txt'" \
            src/selection.slim
done
done
done
done
done
