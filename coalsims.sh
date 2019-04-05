#!/usr/bin/env bash

mkdir -p data/coalsims

for split in `seq 300000 25000 800000`; do
for rep in {1..10}; do
  N=archsplit${split}_afrsplit250000_rep${rep}
  qsub -V -b yes -cwd -j y -l virtual_free=1G,h_vmem=1G -o tmp/sge/${N}.out -N $N \
    ./src/coalsim.py --split-arch ${split} --split-afr 250000 \
      --arch-ages 130000 --ui-age 45000 --neur 5 --nafr 5 --nasn 5 \
      --format snp --output data/coalsims/${N}.tsv
done
done
