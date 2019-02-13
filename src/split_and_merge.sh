#!/usr/bin/env bash

sample_name=$1
bam=$2
index_file=`realpath $3`

mkdir -p tmp/${sample_name}
cd tmp/${sample_name}

/home/mmeyer/perlscripts/solexa/filework/splitBAM.pl -byfile $index_file -overwrite_rg $bam > index_stats.txt

samtools merge ${sample_name}.bam *.bam
