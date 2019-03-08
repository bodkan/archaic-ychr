SHELL := /bin/bash

# directories
data_dir := data
bam_dir := $(data_dir)/bam
fasta_dir := $(data_dir)/fasta
vcf_dir := $(data_dir)/vcf
coord_dir := $(data_dir)/coord
sim_dir := $(data_dir)/sim
tmp_dir := tmp
fig_dir := fig
src_dir := src
dirs := $(data_dir) $(bam_dir) $(vcf_dir) $(fasta_dir) $(coord_dir) $(fig_dir) $(tmp_dir) $(sim_dir) $(tmp_dir)/sge

# BAM files
sgdp_bams := S_BedouinB-1.bam S_Turkish-1.bam S_French-1.bam S_Burmese-1.bam S_Thai-1.bam S_Finnish-2.bam S_Sardinian-1.bam S_Han-2.bam S_Dai-2.bam S_Punjabi-1.bam S_Saami-2.bam S_Papuan-2.bam S_Karitiana-1.bam S_Dinka-1.bam S_Mbuti-1.bam S_Yoruba-2.bam S_Gambian-1.bam S_Mandenka-1.bam S_Ju_hoan_North-1.bam
published_bams := $(sgdp_bams) ustishim.bam a00.bam a00_1.bam a00_2.bam kk1.bam mota.bam bichon.bam loschbour.bam
full_bams := $(addprefix $(bam_dir)/, $(addprefix full_, spy1.bam mez2.bam comb_neand.bam denisova8.bam $(published_bams)))
lippold_bams := $(addprefix $(bam_dir)/, $(addprefix lippold_, elsidron2.bam denisova8.bam $(published_bams)))
exome_bams := $(addprefix $(bam_dir)/, $(addprefix exome_, elsidron1.bam $(published_bams)))

test_bams := $(bam_dir)/control_vindija.bam $(bam_dir)/control_stuttgart.bam

# VCF files
published_vcfs := $(subst .bam,.vcf.gz, $(published_bams))
full_vcfs := $(addprefix $(vcf_dir)/, $(addprefix full_, spy1.vcf.gz mez2.vcf.gz comb_neand.vcf.gz denisova8.vcf.gz $(published_vcfs)))
lippold_vcfs := $(addprefix $(vcf_dir)/, $(addprefix lippold_, elsidron2.vcf.gz denisova8.vcf.gz $(published_vcfs)))
exome_vcfs := $(addprefix $(vcf_dir)/, $(addprefix exome_, elsidron1.vcf.gz $(published_vcfs)))

full_vcf := $(vcf_dir)/merged_full.vcf.gz
full_tv_vcf := $(vcf_dir)/merged_full_tv.vcf.gz
lippold_vcf := $(vcf_dir)/merged_lippold.vcf.gz
exome_vcf := $(vcf_dir)/merged_exome.vcf.gz

test_vcfs := $(vcf_dir)/test_cov.vcf.gz $(vcf_dir)/test_gt.vcf.gz


# FASTA files
fastas := chimp_nea_bteam.fa sidron_bteam.fa bteam.fa
exome_fastas := $(addprefix $(fasta_dir)/exome_,$(fastas))

# scripts
bam_sample := /mnt/expressions/mp/bam-sample/bam-sample
run_nb := $(src_dir)/run_nb.sh
split_and_merge := $(src_dir)/split_and_merge.sh
split_bam := /r1/people/mmeyer/perlscripts/solexa/filework/splitBAM.pl
analyze_bam := /home/mmeyer/perlscripts/solexa/analysis/analyzeBAM.pl

# coordinates
full_bed := $(coord_dir)/capture_full.bed
lippold_bed := $(coord_dir)/capture_lippold.bed
exome_bed := $(coord_dir)/capture_exome.bed
annot_bed := $(coord_dir)/cds.bed $(coord_dir)/phastcons.bed $(coord_dir)/genes.bed $(coord_dir)/pseudogenes.bed
full_sites := $(coord_dir)/capture_full.pos
lippold_sites := $(coord_dir)/capture_lippold.pos
exome_sites := $(coord_dir)/capture_exome.pos

ref_genome := /mnt/solexa/Genomes/hg19_evan/whole_genome.fa

map_filter := /mnt/454/HCNDCAM/Hengs_Alignability_Filter/hs37m_filt35_50.bed.gz


.PHONY: default init bam vcf fasta diagnostics clean



default:
	@echo -e "Usage:"
	@echo -e "\tmake init         -- create all necessary directories"
	@echo -e "\tmake bam          -- process and filter BAM files"
	@echo -e "\tmake vcf          -- run consensus-based genotyping"
	@echo -e "\tmake fasta        -- generate FASTA alignments from VCF files"
	@echo -e "\tmake diagnostics  -- generate diagnostic plots on BAMs"
	@echo -e "\tmake clean        -- delete all generated output file"

init: $(dirs) $(full_bed) $(lippold_bed) $(exome_bed) $(full_sites) $(lippold_sites) $(exome_sites)

bam: $(dirs) $(full_bams) $(lippold_bams) $(exome_bams) $(test_bams)

vcf: $(full_vcf) $(lippold_vcf) $(exome_vcf) $(full_tv_vcf) $(test_vcfs)

fasta: $(dirs) $(full_fastas) $(lippold_fastas) $(exome_fastas)

diagnostics:
	bams=`cd $(bam_dir); ls *French* *a00.bam *elsidron* *ustishim* *kk1* *loschbour* *bichon* *mota* *denisova* *spy* *mez* | grep 'bam$$' | xargs realpath`; \
	mkdir -p $(fig_dir)/damage; \
	cd $(fig_dir)/damage; \
	for b in $$bams; do \
	    qsub -V -cwd -j y -l virtual_free=500M,h_vmem=500M,class=cow /home/mmeyer/perlscripts/solexa/analysis/substitution_patterns.pl $$b; \
	done



#
# BAM processing
#
$(bam_dir)/full_%.bam: $(tmp_dir)/%.bam
	bedtools intersect -a $< -b $(coord_dir)/capture_full.bed > $@
	samtools index $@

$(bam_dir)/lippold_%.bam: $(tmp_dir)/%.bam
	bedtools intersect -a $< -b $(coord_dir)/capture_lippold.bed > $@
	samtools index $@

$(bam_dir)/exome_%.bam: $(tmp_dir)/%.bam
	bedtools intersect -a $< -b $(coord_dir)/capture_exome.bed > $@
	samtools index $@

$(addprefix $(tmp_dir)/, $(sgdp_bams)):
	samtools view -hb --min-tlen 35 -q 25 /mnt/genotyping/sk_pipelines/datasets/Mallick2016_SGDP_Ychromosome/$(basename $(notdir $@)).Y.bam -o $@

# A00 Y
$(tmp_dir)/a00.bam: $(tmp_dir)/a00_1.bam $(tmp_dir)/a00_2.bam
	samtools merge $@ $^

$(tmp_dir)/a00_1.bam:
	cd $(tmp_dir); wget http://evolbio.ut.ee/chrY/GRC13292545.chrY.bam
	samtools view -hb --min-tlen 35 -q 25 $(tmp_dir)/GRC13292545.chrY.bam -o $@

$(tmp_dir)/a00_2.bam:
	cd $(tmp_dir); wget http://evolbio.ut.ee/chrY/GRC13292546.chrY.bam
	samtools view -hb --min-tlen 35 -q 25 $(tmp_dir)/GRC13292546.chrY.bam -o $@

$(tmp_dir)/elsidron1.bam:
	cd $(tmp_dir); curl -O http://cdna.eva.mpg.de/neandertal/exomes/BAM/Sidron_exome_hg19_1000g_LowQualDeamination.md.bam; \
	    samtools index Sidron_exome_hg19_1000g_LowQualDeamination.md.bam; \
	    samtools view -h -b Sidron_exome_hg19_1000g_LowQualDeamination.md.bam Y -o Sidron_exome_hg19_1000g_LowQualDeamination.md.Y.bam; \
	    $(analyze_bam) -qual 25 -minlength 35 Sidron_exome_hg19_1000g_LowQualDeamination.md.Y.bam
	mv $(tmp_dir)/Sidron_exome_hg19_1000g_LowQualDeamination.md.Y.uniq.L35MQ25.bam $@
	samtools index $@

$(tmp_dir)/elsidron2.bam: $(tmp_dir)/elsidron_run1/elsidron_run1.bam $(tmp_dir)/elsidron_run2/elsidron_run2.bam
	samtools merge $(tmp_dir)/elsidron_both_runs.bam $(tmp_dir)/elsidron_run1/elsidron_run1.bam $(tmp_dir)/elsidron_run2/elsidron_run2.bam
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 elsidron_both_runs.bam
	mv $(tmp_dir)/elsidron_both_runs.uniq.L35MQ25.bam $@
	samtools index $@

$(tmp_dir)/elsidron_run1/elsidron_run1.bam:
	$(split_and_merge) elsidron_run1 /mnt/ngs_data/130917_SN7001204_0228_BH06Y0ADXX_R_PEdi_A3207_A3208/Ibis/BWA/s_2-hg19_evan.bam input/A2970_A3206_A3208.txt
$(tmp_dir)/elsidron_run2/elsidron_run2.bam:
	$(split_and_merge) elsidron_run2 /mnt/ngs_data/131129_SN7001204_0235_BH72E4ADXX_R_PEdi_A3601_A3605/Bustard/BWA/s_2_sequence_ancient_hg19_evan.bam input/A2970_A3206_A3208.txt

$(tmp_dir)/ustishim.bam:
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 /mnt/454/Vindija/high_cov/final_bam/Ust_Ishim/chrY.bam
	mv $(tmp_dir)/chrY.uniq.L35MQ25.bam $@
	samtools index $@

$(tmp_dir)/kk1.bam:
	samtools view -h -b /mnt/expressions/mp/Archive/y-selection/tmp/KK1_sort_rmdup_merge_IR_q30_mapDamage.bam Y -o $(tmp_dir)/KK1.Y.bam
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 KK1.Y.bam
	samtools calmd $(tmp_dir)/KK1.Y.uniq.L35MQ25.bam $(ref_genome) --output-fmt BAM > $@
	samtools index $@

$(tmp_dir)/bichon.bam:
	samtools view -h -b /mnt/expressions/mp/Archive/y-selection/tmp/Bichon.sort.rmdup.IR.q30.mapDamage.bam Y -o $(tmp_dir)/Bichon.Y.bam
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 Bichon.Y.bam
	samtools calmd $(tmp_dir)/Bichon.Y.uniq.L35MQ25.bam $(ref_genome) --output-fmt BAM > $@
	samtools index $@

$(tmp_dir)/mota.bam:
	samtools view -h -b /mnt/expressions/mp/Archive/y-selection/tmp/GB20_sort_merge_dedup_l30_IR_q30_mapDamage.bam Y -o $(tmp_dir)/Mota.Y.bam
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 Mota.Y.bam
	samtools calmd $(tmp_dir)/Mota.Y.uniq.L35MQ25.bam $(ref_genome) --output-fmt BAM > $@
	samtools index $@

$(tmp_dir)/loschbour.bam:
	samtools view -h -b /mnt/expressions/mp/Archive/y-selection/tmp/Loschbour.hg19_1000g.bam Y -o $(tmp_dir)/Loschbour.Y.bam
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 Loschbour.Y.bam
	mv $(tmp_dir)/Loschbour.Y.uniq.L35MQ25.bam $@
	samtools index $@

$(tmp_dir)/denisova8.bam:
	$(split_and_merge) denisova8 /mnt/ngs_data/180503_D00829_0138_BCC49NANXX_PEdi_SN_EE_BN_MG/Bustard/BWA/proc1/s_7_sequence_ancient_hg19_evan.bam input/20190207_Ychromosome_Denisova8.txt
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 denisova8/denisova8.bam
	mv $(tmp_dir)/denisova8.uniq.L35MQ25.bam $@
	samtools index $@

$(tmp_dir)/spy1.bam:
	$(split_and_merge) spy1 /mnt/ngs_data/170825_D00829_0064_AHNVL5BCXY_R_PEdi_F5281_F5282/Bustard/BWA/proc1/s_2_sequence_ancient_hg19_evan.bam input/20190207_Ychromosome_Spy1.txt
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 spy1/spy1.bam
	mv $(tmp_dir)/spy1.uniq.L35MQ25.bam $@
	samtools index $@

$(tmp_dir)/mez2.bam:
	$(split_and_merge) mez2 /mnt/ngs_data/170825_D00829_0064_AHNVL5BCXY_R_PEdi_F5281_F5282/Bustard/BWA/proc1/s_2_sequence_ancient_hg19_evan.bam input/20190207_Ychromosome_Mez2.txt
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 mez2/mez2.bam
	mv $(tmp_dir)/mez2.uniq.L35MQ25.bam $@
	samtools index $@

$(tmp_dir)/comb_neand.bam: $(tmp_dir)/mez2.bam $(tmp_dir)/spy1.bam
	samtools merge $@ $^

$(bam_dir)/control_vindija.bam:
	samtools view -h -b /mnt/sequencedb/AncientGenomes/Unpublished/Vi33.19/final_bam/IndelRealign/Vi33.19.chrY.indel_realn.bam Y -o $(tmp_dir)/control_vindija.Y.bam
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 control_vindija.Y.bam
	mv $(tmp_dir)/control_vindija.Y.uniq.L35MQ25.bam $@
	samtools index $@

$(bam_dir)/control_stuttgart.bam:
	samtools view -h -b /mnt/454/Vindija/high_cov/final_bam/LBK/chrY-reali.bam Y -o $(tmp_dir)/control_stuttgart.Y.bam
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 control_stuttgart.Y.bam
	mv $(tmp_dir)/control_stuttgart.Y.uniq.L35MQ25.bam $@
	samtools index $@


#
# VCF processing
#

$(vcf_dir)/%_tv.vcf.gz: $(vcf_dir)/%.vcf.gz
	zless $< | awk -F '\t' '!(($$4 == "A" && $$5 == "G") || ($$4 == "G" && $$5 == "A") || ($$4 == "C" && $$5 == "T") || ($$4 == "T" && $$5 == "C")) || $$0 ~ /^#/' | bgzip > $@
	tabix $@

$(vcf_dir)/merged_full.vcf.gz: $(vcf_dir)/full_chimp.vcf.gz $(full_vcfs)
	bcftools merge $^ | bcftools annotate -x INFO | bcftools view -M 2 -Oz -o $@.all
	bedtools intersect -header -a $@.all -b $(coord_dir)/capture_full.bed | bgzip -c > $@; rm $@.all
	tabix $@

$(vcf_dir)/merged_lippold.vcf.gz: $(vcf_dir)/lippold_chimp.vcf.gz $(lippold_vcfs)
	bcftools merge $^ | bcftools annotate -x INFO | bcftools view -M 2 -Oz -o $@.all
	bedtools intersect -header -a $@.all -b $(coord_dir)/capture_lippold.bed | bgzip -c > $@; rm $@.all
	tabix $@

$(vcf_dir)/merged_exome.vcf.gz: $(vcf_dir)/exome_chimp.vcf.gz $(exome_vcfs)
	bcftools merge $^ |  bcftools annotate -x INFO | bcftools view -M 2 -Oz -o $@.all
	bedtools intersect -header -a $@.all -b $(coord_dir)/capture_exome.bed | bgzip -c > $@; rm $@.all
	tabix $@

# generate genotypes from the Chimp reference genome
$(vcf_dir)/%_chimp.vcf.gz: $(coord_dir)/capture_%.pos
	$(src_dir)/chimp_vcf.sh $< $(basename $@)
	bgzip $(basename $@)
	tabix $@

# genotype samples by consensus calling
$(vcf_dir)/%.vcf.gz: $(bam_dir)/%.bam
	$(bam_sample) --bam $< --ref $(ref_genome) --format vcf \
	    --strategy consensus --mincov 4 --minbq 20 --minmq 25 \
	    --sample-name $(shell echo $(basename $(notdir $<)) | sed 's/^[a-z]*_//') --output $(basename $(basename $@))
	bgzip $(basename $@)
	tabix $@

#
# testing the coverage cutoff
#
test_cov_vcfs := $(addsuffix .vcf.gz, $(addprefix $(tmp_dir)/, $(addprefix test_cov_,$(shell seq 1 9))))

$(vcf_dir)/test_cov.vcf.gz: $(vcf_dir)/full_a00.vcf.gz $(test_cov_vcfs)
	bcftools merge $^ | bcftools view -M 2 -Oz -o $@
	tabix $@

$(tmp_dir)/test_cov_%.vcf.gz: $(bam_dir)/full_denisova8.bam
	depth=`echo $@ | sed 's/.*\([0-9]\).*/\1/'`; \
	$(bam_sample) --bam $< --ref $(ref_genome) --format vcf \
	    --strategy consensus --mincov $$depth --minbq 20 --minmq 25 \
	    --sample-name $(shell echo $(basename $(notdir $<)) | sed 's/^[a-z]*_//')_dp$${depth} --output $(basename $(basename $@))
	bgzip $(basename $@)
	tabix $@

#
# testing A00 VCF file for comparing bam-sample and bcftools calls
#
$(vcf_dir)/test_gt.vcf.gz: $(vcf_dir)/full_a00.vcf.gz $(vcf_dir)/full_denisova8.vcf.gz $(vcf_dir)/full_ustishim.vcf.gz
	samtools mpileup -B -t DP -Q 20 -q 25 -u -f $(ref_genome) $(bam_dir)/full_a00.bam \
		| bcftools call --ploidy 1 -m -V indels -Oz -o $(tmp_dir)/bcftools_a00.vcf.gz
	tabix $(tmp_dir)/bcftools_a00.vcf.gz
	samtools mpileup -B -t DP -Q 20 -q 25 -u -f $(ref_genome) $(bam_dir)/full_denisova8.bam \
		| bcftools call --ploidy 1 -m -V indels -Oz -o $(tmp_dir)/bcftools_denisova8.vcf.gz
	tabix $(tmp_dir)/bcftools_denisova8.vcf.gz
	samtools mpileup -B -t DP -Q 20 -q 30 -u -f $(ref_genome) $(bam_dir)/full_ustishim.bam \
		| bcftools call --ploidy 1 -m -V indels -Oz -o $(tmp_dir)/bcftools_ustishim.vcf.gz
	tabix $(tmp_dir)/bcftools_ustishim.vcf.gz
	bcftools merge $(vcf_dir)/full_a00.vcf.gz $(tmp_dir)/bcftools_a00.vcf.gz $(vcf_dir)/full_denisova8.vcf.gz $(tmp_dir)/bcftools_denisova8.vcf.gz $(vcf_dir)/full_ustishim.vcf.gz $(tmp_dir)/bcftools_ustishim.vcf.gz \
	    | bcftools view -v snps \
	    | bcftools annotate -x INFO,FORMAT/PL \
	    | bcftools reheader -s <(echo -e "consensus_a00\nbcftools_a00\nconsensus_denisova8\nbcftools_denisova8\nconsensus_ustishim\nbcftools_ustishim\n"| cat) \
	    | bgzip -c > $@
	tabix $@



#
# FASTA alignments for BEAST analyses
#
sample_ids := ElSidron A00 $(bteam_names)

$(fasta_dir)/%_chimp_nea_bteam.fa: $(vcf_dir)/merged_%.vcf.gz
	python $(src_dir)/vcf_to_fasta.py --vcf-file $< --fasta-file $@ --chrom Y --sample-names Chimp Mez2 Spy $(sample_ids)

$(fasta_dir)/%_sidron_bteam.fa: $(vcf_dir)/merged_%.vcf.gz
	python $(src_dir)/vcf_to_fasta.py --vcf-file $< --fasta-file $@ --chrom Y --sample-names $(sample_ids)

$(fasta_dir)/%_bteam.fa: $(vcf_dir)/merged_%.vcf.gz
	python $(src_dir)/vcf_to_fasta.py --vcf-file $< --fasta-file $@ --chrom Y --sample-names A00 $(bteam_names)



#
# coordinate files
#

# Y chromosome capture regions from Lippold et al. (~570 kb)
# /mnt/genotyping/sendru/basti_design.bed
$(lippold_bed):
	# cp input/basti_design.bed > $@
	bedtools intersect -a input/basti_design.bed -b $(map_filter) > $@.tmp
	bedtools sort -i $@.tmp > $@; rm $@.tmp


# Y chromosome capture regions designed by Qiaomei (~6Mb)
$(full_bed):
	# perl -lane 'print $$F[0] . "\t" . $$F[1] . "\t" . $$F[2]' input/Y.filt35_50_SRepeat_100.bed > $@
	bedtools intersect \
	    -a <(perl -lane 'print $$F[0] . "\t" . $$F[1] . "\t" . $$F[2]' input/Y.filt35_50_SRepeat_100.bed) \
	    -b $(map_filter) \
	    > $@

# Y chromosome exome capture regions
$(exome_bed):
	cd $(tmp_dir); wget http://www.cell.com/cms/attachment/2052899616/2060015784/mmc2.zip; unzip mmc2.zip
	cp $(tmp_dir)/ajhg2064mmc2_V1.txt $@
	# bedtools intersect -a $(tmp_dir)/ajhg2064mmc2_V1.txt -b $(map_filter) > $@

# functional annotation coordinates
$(annot_bed):
	$(run_nb) notebooks/annotated_regions.ipynb

# sites within BED files
$(coord_dir)/capture_%.pos: $(coord_dir)/capture_%.bed
	$(src_dir)/sites_in_bed.py --bed $< --format pos --output $@



$(dirs):
	mkdir -p $@

clean:
	rm -rf $(dirs)
