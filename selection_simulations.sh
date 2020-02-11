#!/usr/bin/env bash

# simulations exploring the admixture time vs amount of sequence
direction="modern2neand"
mkdir -p data/sim/${direction}

for gene_total in `seq 100000 100000 2000000`; do
for admix_rate in 0.05; do
for admix_time in `seq 175000 25000 450000`; do
for rep in `seq 1 100`; do
    N="seq${gene_total}_time${admix_time}_rate${admix_rate}_rep${rep}"
    qsub -V -b yes -cwd -j y -l virtual_free=5G,h_vmem=5G -o tmp/sge/${N}.out -N $N \
        slim \
            -d mut_rate=1.85e-08 \
            -d gene_total=${gene_total} \
            -d "direction='${direction}'" \
            -d admix_time=${admix_time} \
            -d admix_rate=${admix_rate} \
            -d dump_at=$((${admix_time} + 25)) \
            -d sample_for=100000 \
            -d "output='data/sim/${direction}/${N}'" \
            src/selection.slim
done
done
done
done
