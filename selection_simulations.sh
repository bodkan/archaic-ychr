#!/usr/bin/env bash

for type in A X Y; do
for rep in `seq 1 100`; do
    N="archaic_chr${type}_sel_rep${rep}"
    qsub -V -b yes -cwd -j y -l virtual_free=2G,h_vmem=2G,class=cow -o tmp/sge/${N}.out -N $N \
        slim \
            -d "chrom_type='${type}'" \
            -d mut_rate=1e-8 \
            -d chrom_length=60000000 \
            -d gene_total=3000000 \
            -d gene_count=60 \
            -d admix_time=55000 \
            -d archaic_rate=0.1 \
            -d modern_rate=0 \
            -d "output='data/sim/${N}.txt'" \
            src/admixture.slim
done
done

# slim -d 'chrom_type="A"' \
#      -d mut_rate=1e-8 \
#      -d chrom_length=60000000 \
#      -d gene_total=3000000 \
#      -d gene_count=60 \
#      -d admix_time=55000 \
#      -d archaic_rate=0.1 \
#      -d modern_rate=0 \
#      -d 'output="asd.txt"' \
#      src/admixture.slim

