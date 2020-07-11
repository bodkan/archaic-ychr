#!/usr/bin/env bash

mkdir /mnt/scratch/mp/zenodo
cd /mnt/scratch/mp/zenodo

cp /mnt/454/Carbon_beast_QM/array_2015_0729/array_order/Y.filt35_50_SRepeat_100.bed .
cp /mnt/454/Carbon_beast_QM/array_2015_0729/Y.filt35_50_SRepeat_100_3tiling.probes .

cp /mnt/454/array_David/basti_chrY_region.bed .
cp /mnt/454/array_David/basti_chrY_region.probes .

mv Y.filt35_50_SRepeat_100.bed target_Y_whole.bed
mv Y.filt35_50_SRepeat_100_3tiling.probes target_Y_whole.probes

mv basti_chrY_region.bed target_Y_560kb.bed
mv basti_chrY_region.probes target_Y_560kb.probes

tar cvf capture_designs.tar.gz *.bed *.probes

