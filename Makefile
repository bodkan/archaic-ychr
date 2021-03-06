SHELL := /bin/bash

# directories
data_dir := data
bam_dir := $(data_dir)/bam
pileup_dir := $(data_dir)/pileup
fasta_dir := $(data_dir)/fasta
vcf_dir := $(data_dir)/vcf
coord_dir := $(data_dir)/coord
sim_dir := $(data_dir)/sim
rds_dir := $(data_dir)/rds
test_dir := test
tmp_dir := tmp
fig_dir := figures
src_dir := src
dirs := $(data_dir) $(bam_dir) $(pileup_dir) $(vcf_dir) $(fasta_dir) $(test_dir) $(coord_dir) $(fig_dir) $(tmp_dir) $(sim_dir) $(rds_dir) $(tmp_dir)/sge

# BAM files
mez2_subsamples := $(addprefix mez2_dp, $(shell seq 1 10))
elsidron2_subsamples := $(addprefix elsidron2_dp, $(shell seq 1 7))
sgdp_bams := S_BedouinB-1.bam S_Turkish-1.bam S_French-1.bam S_Burmese-1.bam S_Thai-1.bam S_Finnish-2.bam S_Sardinian-1.bam S_Han-2.bam S_Dai-2.bam S_Punjabi-1.bam S_Saami-2.bam S_Papuan-2.bam S_Karitiana-1.bam S_Ju_hoan_North-1.bam S_Dinka-1.bam S_Mbuti-1.bam S_Yoruba-2.bam S_Gambian-1.bam S_Mandenka-1.bam
modern_bams := a00.bam a00_1.bam a00_2.bam $(sgdp_bams)
full_bams := $(addprefix $(bam_dir)/, $(addprefix full_, ustishim.bam shotgun_spy1.bam shotgun_mez2.bam spy1.bam mez2.bam den8.bam den4.bam den_merged.bam spy1_merged.bam $(addsuffix .bam, $(mez2_subsamples)) $(modern_bams)))
lippold_bams := $(addprefix $(bam_dir)/, $(addprefix lippold_, ustishim.bam spy1.bam mez2.bam elsidron2.bam den8.bam den4.bam spy1_merged.bam $(addsuffix .bam, $(elsidron2_subsamples)) $(modern_bams)))
exome_bams := $(addprefix $(bam_dir)/, $(addprefix exome_, ustishim.bam spy1.bam mez2.bam elsidron1.bam den8.bam den4.bam spy1_merged.bam $(modern_bams)))

test_bams := $(bam_dir)/control_vindija.bam $(bam_dir)/control_stuttgart.bam

# pileup files
pileups := $(addprefix $(pileup_dir)/, $(addprefix full_, spy1.txt.gz mez2.txt.gz den8.txt.gz den4.txt.gz S_French-1.txt.gz a00.txt.gz S_Saami-2.txt.gz S_Han-2.txt.gz)) $(pileup_dir)/lippold_elsidron2.txt.gz

# VCF files
modern_vcfs := $(subst .bam,.vcf.gz, $(modern_bams))
full_arch_vcfs    := $(addprefix $(vcf_dir)/, $(addprefix full_, ustishim.vcf.gz shotgun_spy1.vcf.gz shotgun_mez2.vcf.gz spy1.vcf.gz mez2.vcf.gz den8.vcf.gz den4.vcf.gz den_merged.vcf.gz spy1_merged.vcf.gz mez2_snpad.vcf.gz spy1_snpad.vcf.gz den4_snpad.vcf.gz den8_snpad.vcf.gz $(addsuffix .vcf.gz, $(mez2_subsamples))))
lippold_arch_vcfs := $(addprefix $(vcf_dir)/, $(addprefix lippold_, ustishim.vcf.gz spy1.vcf.gz mez2.vcf.gz elsidron2.vcf.gz den8.vcf.gz den4.vcf.gz spy1_merged.vcf.gz $(addsuffix .vcf.gz, $(elsidron2_subsamples))))
exome_arch_vcfs   := $(addprefix $(vcf_dir)/, $(addprefix exome_, ustishim.vcf.gz spy1.vcf.gz mez2.vcf.gz elsidron1.vcf.gz den8.vcf.gz den4.vcf.gz spy1_merged.vcf.gz))
full_modern_vcfs     := $(addprefix $(vcf_dir)/, $(addprefix full_, $(modern_vcfs)))
lippold_modern_vcfs  := $(addprefix $(vcf_dir)/, $(addprefix lippold_, $(modern_vcfs)))
exome_modern_vcfs    := $(addprefix $(vcf_dir)/, $(addprefix exome_, $(modern_vcfs)))

full_vcf := $(vcf_dir)/full_modern.vcf.gz
lippold_vcf := $(vcf_dir)/lippold_modern.vcf.gz
exome_vcf := $(vcf_dir)/exome_modern.vcf.gz

test_vcfs := $(foreach sample, spy1 den4 den8 mez2 a00, $(test_dir)/genotyping_$(sample).vcf.gz)

# FASTA files
fastas := $(addprefix $(fasta_dir)/,full_merged_nodmg.fa full_merged_allsnps.fa modern_all_full_merged.fa modern_var_full_merged.fa nochimp_full_merged_nodmg.fa nochimp_full_merged_allsnps.fa highcov_full_merged_allsnps.fa highcov_full_merged_nodmg.fa)

# scripts
bam_caller := $(src_dir)/bam-caller.py
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
	@echo -e "\tmake pileup       -- generate pileup data from BAM files"
	@echo -e "\tmake vcf          -- run consensus-based genotyping"
	@echo -e "\tmake fasta        -- generate FASTA alignments from VCF files"
	@echo -e "\tmake diagnostics  -- generate diagnostic plots on BAMs"
	@echo -e "\tmake simulations  -- run selection simulations using SLiM"
	@echo -e "\tmake clean        -- delete all generated output file"

init: $(dirs) $(full_bed) $(lippold_bed) $(exome_bed) $(full_sites) $(lippold_sites) $(exome_sites)

bam: $(full_bams) $(lippold_bams) $(exome_bams) $(test_bams)

pileup: $(pileups)

vcf: $(full_arch_vcfs) $(lippold_arch_vcfs) $(exome_arch_vcfs) $(full_vcf) $(lippold_vcf) $(exome_vcf) $(test_vcfs)

fasta: $(fastas)

diagnostics:
	bams=`cd $(bam_dir); ls exome_elsidron1.bam lippold_elsidron2.bam full_den4.bam full_den8.bam full_spy1.bam full_mez2.bam full_shotgun*.bam | grep 'bam$$' | xargs realpath`; \
	mkdir -p $(data_dir)/damage; \
	cd $(data_dir)/damage; \
	for b in $$bams; do \
		/home/mmeyer/perlscripts/solexa/analysis/substitution_patterns.pl $$b & \
	done

simulations:
	./src/selection_simulations.sh



######################################################################
# BAM processing
######################################################################

$(bam_dir)/full_mez2_dp%.bam: $(bam_dir)/full_mez2.bam $(coord_dir)/capture_full.bed
	target_coverage=`basename $@ .bam | sed 's/full_mez2_dp//'`; \
	input_coverage=`bedtools coverage -a $(coord_dir)/capture_full.bed -b $< -d | awk '{sum+=$$5} END { print sum/NR}'`; \
	samtools view -b -s `echo $$target_coverage $$input_coverage | awk '{print $$1/$$2}'` $< > $@
	samtools index $@

$(bam_dir)/lippold_elsidron2_dp%.bam: $(bam_dir)/lippold_elsidron2.bam $(coord_dir)/capture_lippold.bed
	target_coverage=`basename $@ .bam | sed 's/lippold_elsidron2_dp//'`; \
	input_coverage=`bedtools coverage -a $(coord_dir)/capture_lippold.bed -b $< -d | awk '{sum+=$$5} END { print sum/NR}'`; \
	samtools view -b -s `echo $$target_coverage $$input_coverage | awk '{print $$1/$$2}'` $< > $@
	samtools index $@

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
	samtools index $@

$(tmp_dir)/a00_1.bam:
	cd $(tmp_dir); wget http://evolbio.ut.ee/chrY/GRC13292545.chrY.bam
	samtools view -hb --min-tlen 35 -q 25 $(tmp_dir)/GRC13292545.chrY.bam -o $@
	samtools index $@

$(tmp_dir)/a00_2.bam:
	cd $(tmp_dir); wget http://evolbio.ut.ee/chrY/GRC13292546.chrY.bam
	samtools view -hb --min-tlen 35 -q 25 $(tmp_dir)/GRC13292546.chrY.bam -o $@
	samtools index $@

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

$(tmp_dir)/den8.bam:
	$(split_and_merge) den8 /mnt/ngs_data/180503_D00829_0138_BCC49NANXX_PEdi_SN_EE_BN_MG/Bustard/BWA/proc1/s_7_sequence_ancient_hg19_evan.bam input/20190207_Ychromosome_Denisova8.txt
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 den8/den8.bam
	mv $(tmp_dir)/den8.uniq.L35MQ25.bam $@
	samtools index $@

$(tmp_dir)/den_merged.bam: $(tmp_dir)/den4.bam $(tmp_dir)/den8.bam
	samtools merge $@ $^
	samtools index $@

$(tmp_dir)/spy1_merged.bam: $(tmp_dir)/spy1.bam $(tmp_dir)/shotgun_spy1.bam
	samtools merge $@ $^
	samtools index $@

$(tmp_dir)/spy1.bam:
	$(split_and_merge) spy1 /mnt/ngs_data/170825_D00829_0064_AHNVL5BCXY_R_PEdi_F5281_F5282/Bustard/BWA/proc1/s_2_sequence_ancient_hg19_evan.bam input/20190207_Ychromosome_Spy1.txt
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 spy1/spy1.bam
	mv $(tmp_dir)/spy1.uniq.L35MQ25.bam $@
	samtools index $@

# Mezmaiskaya first capture
$(tmp_dir)/mez2_capture1.bam:
	$(split_and_merge) mez2_capture1 /mnt/ngs_data/170825_D00829_0064_AHNVL5BCXY_R_PEdi_F5281_F5282/Bustard/BWA/proc1/s_2_sequence_ancient_hg19_evan.bam input/20190207_Ychromosome_Mez2.txt
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 mez2_capture1/mez2_capture1.bam
	mv $(tmp_dir)/mez2_capture1.uniq.L35MQ25.bam $@
	samtools index $@

# Mezmaiskaya second capture - lane 1
$(tmp_dir)/mez2_capture2_lane1.bam:
	$(split_and_merge) mez2_capture2_lane1 /mnt/scratch/janet/illumina/190403_D00829_0243_BHY353BCX2_R_PEdi_N3550_N3551_reprocess/Bustard/BWA/proc2/s_1_sequence_ancient_hg19_evan.bam input/20190416_Mez2_s1.txt
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 mez2_capture2_lane1/mez2_capture2_lane1.bam
	mv $(tmp_dir)/mez2_capture2_lane1.uniq.L35MQ25.bam $@
	samtools index $@

# Mezmaiskaya second capture - lane 2
$(tmp_dir)/mez2_capture2_lane2.bam:
	$(split_and_merge) mez2_capture2_lane2 /mnt/scratch/janet/illumina/190403_D00829_0243_BHY353BCX2_R_PEdi_N3550_N3551_reprocess/Bustard/BWA/proc2/s_2_sequence_ancient_hg19_evan.bam input/20190416_Mez2_s2.txt
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 mez2_capture2_lane2/mez2_capture2_lane2.bam
	mv $(tmp_dir)/mez2_capture2_lane2.uniq.L35MQ25.bam $@
	samtools index $@

# merge of all Mez2 samples
$(tmp_dir)/mez2.bam: $(tmp_dir)/mez2_capture1.bam $(tmp_dir)/mez2_capture2_lane1.bam $(tmp_dir)/mez2_capture2_lane2.bam
	samtools merge $@ $^

# Denisova 4 - lane 1
$(tmp_dir)/den4_capture2_lane1.bam:
	$(split_and_merge) den4_capture2_lane1 /mnt/scratch/janet/illumina/190403_D00829_0243_BHY353BCX2_R_PEdi_N3550_N3551_reprocess/Bustard/BWA/proc2/s_1_sequence_ancient_hg19_evan.bam input/20190416_Den4_s1.txt
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 den4_capture2_lane1/den4_capture2_lane1.bam
	mv $(tmp_dir)/den4_capture2_lane1.uniq.L35MQ25.bam $@
	samtools index $@

# Denisova 4 - lane 2
$(tmp_dir)/den4_capture2_lane2.bam:
	$(split_and_merge) den4_capture2_lane2 /mnt/scratch/janet/illumina/190403_D00829_0243_BHY353BCX2_R_PEdi_N3550_N3551_reprocess/Bustard/BWA/proc2/s_2_sequence_ancient_hg19_evan.bam input/20190416_Den4_s2.txt
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 den4_capture2_lane2/den4_capture2_lane2.bam
	mv $(tmp_dir)/den4_capture2_lane2.uniq.L35MQ25.bam $@
	samtools index $@

# merge of both Den4 lanes
$(tmp_dir)/den4.bam: $(tmp_dir)/den4_capture2_lane1.bam $(tmp_dir)/den4_capture2_lane2.bam
	samtools merge $@ $^

$(tmp_dir)/shotgun_mez2.bam:
#	samtools view -h -b /mnt/expressions/mateja/Late_Neandertals/Final_complete_dataset/L35MQ25_per_individual/Mezmaiskaya2_final.L35MQ25.bam Y -o $@
	cp /mnt/expressions/mp/Archive/late_nea_y_bams/shotgun_mez2.bam $@
	samtools index $@

$(tmp_dir)/shotgun_spy1.bam:
#	samtools view -h -b /mnt/expressions/mateja/Late_Neandertals/Final_complete_dataset/L35MQ25_per_individual/Spy_final.L35MQ25.bam Y -o $@
	cp /mnt/expressions/mp/Archive/late_nea_y_bams/shotgun_spy1.bam $@
	samtools index $@

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



######################################################################
# VCF processing
######################################################################

$(vcf_dir)/full_modern.vcf.gz: $(vcf_dir)/full_chimp.vcf.gz $(vcf_dir)/full_ustishim.vcf.gz $(full_modern_vcfs)
	bcftools merge $^ | bcftools annotate -x INFO | bcftools view -M 2 -Oz -o $@.all
	bedtools intersect -header -a $@.all -b $(coord_dir)/capture_full.bed | bgzip -c > $@; rm $@.all
	tabix $@

$(vcf_dir)/lippold_modern.vcf.gz: $(vcf_dir)/lippold_chimp.vcf.gz $(vcf_dir)/full_ustishim.vcf.gz $(lippold_modern_vcfs)
	bcftools merge $^ | bcftools annotate -x INFO | bcftools view -M 2 -Oz -o $@.all
	bedtools intersect -header -a $@.all -b $(coord_dir)/capture_lippold.bed | bgzip -c > $@; rm $@.all
	tabix $@

$(vcf_dir)/exome_modern.vcf.gz: $(vcf_dir)/exome_chimp.vcf.gz $(vcf_dir)/full_ustishim.vcf.gz $(exome_modern_vcfs)
	bcftools merge $^ |  bcftools annotate -x INFO | bcftools view -M 2 -Oz -o $@.all
	bedtools intersect -header -a $@.all -b $(coord_dir)/capture_exome.bed | bgzip -c > $@; rm $@.all
	tabix $@

# full_mez2_snpad.vcf.gz:
# 	cp /mnt/expressions/Janet/YChr/Mez2/chrY.mq25.map50.vcf $(tmp_dir)
# 	bcftools reheader -s <(echo "full_Mez2map50 mez2_snpad") $(tmp_dir)/chrY.mq25.map50.vcf \
# 	    | bgzip -c > $@
# 	tabix $@

$(vcf_dir)/%_snpad.vcf.gz: $(bam_dir)/%.bam
	name="$(shell echo $(basename $(notdir $<)) | sed 's/^[a-z]*_//')_snpad"; \
	$(src_dir)/call_snpad.sh $$name $< $@
	tabix $@

# generate genotypes from the Chimp reference genome
$(vcf_dir)/%_chimp.vcf.gz: $(coord_dir)/capture_%.pos
	$(src_dir)/chimp_vcf.sh $< $(basename $@)
	bgzip $(basename $@)
	tabix $@

# genotype samples by consensus calling
$(vcf_dir)/%.vcf.gz: $(bam_dir)/%.bam
	name="$(shell echo $(basename $(notdir $<)) | sed 's/^[a-z]*_//')"; \
	$(bam_caller) --bam $< \
	    --strategy majority --proportion 0.9 --mincov 1 --minbq 20 --minmq 25 \
	    --sample-name $$name --output $(basename $(basename $@))
	bgzip $(basename $@)
	tabix $@

#
# testing different genotyping procedures
#
.SECONDARY:

$(test_dir)/genotyping_%.vcf.gz: $(test_dir)/%_baq.vcf.gz $(test_dir)/%_nobaq.vcf.gz $(test_dir)/%_consensus.vcf.gz $(test_dir)/%_tolerance.vcf.gz
	bcftools merge $^ \
		| bcftools annotate -x INFO,FORMAT/PL \
		| bgzip -c \
	> $@
	tabix $@
$(test_dir)/%_baq.vcf.gz: $(bam_dir)/full_%.bam
	bcftools mpileup --min-BQ 20 --min-MQ 25 --annotate FORMAT/DP -Ou -f $(ref_genome) $^ \
		| bcftools call --ploidy 1 -m -V indels \
		| bcftools reheader -s <(echo "baq") \
		| bcftools view - -Oz -o $@
	bedtools intersect -header -a $@ -b $(coord_dir)/capture_full.bed | bgzip -c > $@.filt; mv $@.filt $@
	tabix $@
$(test_dir)/%_nobaq.vcf.gz: $(bam_dir)/full_%.bam
	bcftools mpileup --no-BAQ --min-BQ 20 --min-MQ 25 --annotate FORMAT/DP -Ou -f $(ref_genome) $^ \
		| bcftools call --ploidy 1 -m -V indels \
		| bcftools reheader -s <(echo "nobaq") \
		| bcftools view - -Oz -o $@
	bedtools intersect -header -a $@ -b $(coord_dir)/capture_full.bed | bgzip -c > $@.filt; mv $@.filt $@
	tabix $@
$(test_dir)/%_consensus.vcf.gz: $(bam_dir)/full_%.bam
	$(bam_caller) --bam $< \
	    --strategy majority --proportion 1.0 --mincov 1 --minbq 20 --minmq 25 \
	    --sample-name cons --output $(basename $(basename $@))
	bgzip $(basename $@)
	bedtools intersect -header -a $@ -b $(coord_dir)/capture_full.bed | bgzip -c > $@.filt; mv $@.filt $@
	tabix $@
$(test_dir)/%_tolerance.vcf.gz: $(bam_dir)/full_%.bam
	$(bam_caller) --bam $< \
	    --strategy majority --proportion 0.9 --mincov 1 --minbq 20 --minmq 25 \
	    --sample-name tol --output $(basename $(basename $@))
	bgzip $(basename $@)
	bedtools intersect -header -a $@ -b $(coord_dir)/capture_full.bed | bgzip -c > $@.filt; mv $@.filt $@
	tabix $@


#
# FASTA alignments for BEAST analyses
#

$(vcf_dir)/full_merged.vcf.gz: $(foreach sample,den4 den8 spy1 mez2 modern,$(vcf_dir)/full_$(sample).vcf.gz) $(vcf_dir)/lippold_elsidron2.vcf.gz
	mkdir -p $(tmp_dir)/vcf_fasta
	for f in $(vcf_dir)/full_{den4,den8,spy1,mez2,a00}.vcf.gz $(vcf_dir)/full_S_*.vcf.gz $(vcf_dir)/lippold_elsidron2.vcf.gz; do \
		$(src_dir)/filter_vcf.sh $${f}; \
	done
	bcftools merge $(tmp_dir)/vcf_fasta/full_*.vcf.gz $(tmp_dir)/vcf_fasta/lippold_elsidron2.vcf.gz $(vcf_dir)/full_chimp.vcf.gz | bcftools view -M 2 -Oz -o $@.all
	bedtools intersect -header -a $@.all -b $(coord_dir)/capture_full.bed | bgzip -c > $@; rm $@.all
	tabix $@

# FASTAS for nj tree building
$(fasta_dir)/%_allsnps.fa: $(vcf_dir)/%.vcf.gz
	python $(src_dir)/vcf_to_fasta.py --vcf $< --fasta $@ --variable
$(fasta_dir)/%_nodmg.fa: $(vcf_dir)/%.vcf.gz
	python $(src_dir)/vcf_to_fasta.py --vcf $< --fasta $@ --variable --no-damage

# FASTAS without chimp for BEAST analyses
$(fasta_dir)/nochimp_%_allsnps.fa: $(vcf_dir)/%.vcf.gz
	python $(src_dir)/vcf_to_fasta.py --vcf $< --fasta $@ --variable --exclude chimp
$(fasta_dir)/nochimp_%_nodmg.fa: $(vcf_dir)/%.vcf.gz
	python $(src_dir)/vcf_to_fasta.py --vcf $< --fasta $@ --variable --no-damage --exclude chimp

# FASTAS without chimp with only Mez2 and El Sidron 1253 for BEAST analyses
$(fasta_dir)/highcov_%_allsnps.fa: $(vcf_dir)/%.vcf.gz
	python $(src_dir)/vcf_to_fasta.py --vcf $< --fasta $@ --variable --exclude chimp den4 den8 spy1
$(fasta_dir)/highcov_%_nodmg.fa: $(vcf_dir)/%.vcf.gz
	python $(src_dir)/vcf_to_fasta.py --vcf $< --fasta $@ --variable --no-damage --exclude chimp den4 den8 spy1

# FASTA with AMH + chimp for a contamination tree figure
$(fasta_dir)/modern_all_%.fa: $(vcf_dir)/%.vcf.gz
	python $(src_dir)/vcf_to_fasta.py --vcf $< --fasta $@ --exclude mez2 spy1 den4 den8 elsidron2
$(fasta_dir)/modern_var_%.fa: $(vcf_dir)/%.vcf.gz
	python $(src_dir)/vcf_to_fasta.py --vcf $< --fasta $@ --variable --exclude mez2 spy1 den4 den8 elsidron2


#
# generate pileup files for read-based contamination estimates
#
$(pileup_dir)/%.txt.gz: $(bam_dir)/%.bam
	$(bam_caller) --bam $< --strategy pileup \
	    --minbq 20 --minmq 25 \
	    --output $(basename $(basename $@))
	bgzip $(basename $@)




######################################################################
# coordinate files
######################################################################

# Y chromosome capture regions from Lippold et al. (~570 kb)
# /mnt/genotyping/sendru/basti_design.bed
$(lippold_bed):
	# cp input/basti_design.bed > $@
	bedtools intersect -a input/basti_design.bed -b $(map_filter) > $@.tmp
	bedtools sort -i $@.tmp > $@; rm $@.tmp

# Y chromosome capture regions designed by Qiaomei (~6Mb)
# /mnt/454/Carbon_beast_QM/array_2015_0729/array_order/Y.filt35_50_SRepeat_100.bed
$(full_bed):
	# perl -lane 'print $$F[0] . "\t" . $$F[1] . "\t" . $$F[2]' input/Y.filt35_50_SRepeat_100.bed > $@
	bedtools intersect \
	    -a <(perl -lane 'print $$F[0] . "\t" . $$F[1] . "\t" . $$F[2] if $$F[2] < 30000000' input/Y.filt35_50_SRepeat_100.bed) \
	    -b $(map_filter) \
	    > $@

# Y chromosome exome capture regions
$(exome_bed):
	cd $(tmp_dir); wget http://www.cell.com/cms/attachment/2052899616/2060015784/mmc2.zip; unzip mmc2.zip
	cp $(tmp_dir)/ajhg2064mmc2_V1.txt $@
	# bedtools intersect -a $(tmp_dir)/ajhg2064mmc2_V1.txt -b $(map_filter) > $@

# sites within BED files
$(coord_dir)/capture_%.pos: $(coord_dir)/capture_%.bed
	$(src_dir)/sites_in_bed.py --bed $< --format pos --output $@



######################################################################
# other dependencies
######################################################################

$(dirs):
	mkdir -p $@

clean:
	rm -rf $(dirs)
