#!/usr/bin/env bash

mkdir /mnt/scratch/mp/archaic-ychr
cd /mnt/scratch/mp/archaic-ychr

mkdir data

cp -r /mnt/expressions/mp/archaic-ychr/data/{bam,coord,damage,fasta,pileup,sim,vcf} data/

cp /mnt/454/Carbon_beast_QM/array_2015_0729/array_order/Y.filt35_50_SRepeat_100.bed .
cp /mnt/454/Carbon_beast_QM/array_2015_0729/Y.filt35_50_SRepeat_100_3tiling.probes .

tree data/ > data_contents.txt

tar -zcf data.tar.gz data
sha256sum data.tar.gz > data.tar.gz.sha256sum

rm -rf /mnt/scratch/mp/archaic-ychr/data/
