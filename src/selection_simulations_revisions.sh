#!/usr/bin/env bash

# simulations exploring the trade-off between AMH and Neandertal
# effective population size
mkdir -p data/sim/demography

for gene_total in 500000; do # fixing these
for admix_time in 200000; do # two parameters
for human_Ne in `seq 5000 1000 15000`; do
for neand_Ne in `seq 200 200 2000`; do
for admix_rate in `seq 0.01 0.01 0.1`; do
for rep in `seq 1 25`; do
    N="humanNe${human_Ne}_neandNe${neand_Ne}_rate${admix_rate}_rep${rep}"
    qsub -V -b yes -cwd -j y -l virtual_free=5G,h_vmem=5G -o tmp/sge/${N}.out -N $N \
        slim \
            -d mut_rate=1.85e-08 \
            -d human_Ne=${human_Ne} \
            -d neand_Ne=${neand_Ne} \
            -d gene_total=${gene_total} \
            -d "direction='modern2neand'" \
            -d admix_time=${admix_time} \
            -d admix_rate=${admix_rate} \
            -d dump_at=$((${admix_time} + 25)) \
            -d sample_for=100000 \
            -d "output='data/sim/demography/${N}'" \
            src/selection.slim
done
done
done
done
done
done
