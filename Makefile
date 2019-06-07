SHELL := /bin/bash

# directories
data_dir := data
bam_dir := $(data_dir)/bam
pileup_dir := $(data_dir)/pileup
fasta_dir := $(data_dir)/fasta
vcf_dir := $(data_dir)/vcf
coord_dir := $(data_dir)/coord
sim_dir := $(data_dir)/sim
tmp_dir := tmp
fig_dir := fig
src_dir := src
dirs := $(data_dir) $(bam_dir) $(pileup_dir) $(vcf_dir) $(fasta_dir) $(coord_dir) $(fig_dir) $(tmp_dir) $(sim_dir) $(tmp_dir)/sge

# BAM files
sgdp_bams := S_BedouinB-1.bam S_Turkish-1.bam S_French-1.bam S_Burmese-1.bam S_Thai-1.bam S_Finnish-2.bam S_Sardinian-1.bam S_Han-2.bam S_Dai-2.bam S_Punjabi-1.bam S_Saami-2.bam S_Papuan-2.bam S_Karitiana-1.bam S_Dinka-1.bam S_Mbuti-1.bam S_Yoruba-2.bam S_Gambian-1.bam S_Mandenka-1.bam S_Ju_hoan_North-1.bam
published_bams := ustishim.bam a00.bam $(sgdp_bams)
full_bams := $(addprefix $(bam_dir)/, $(addprefix full_, shotgun_spy1.bam shotgun_mez2.bam spy1.bam mez2.bam den8.bam den4.bam $(published_bams)))
lippold_bams := $(addprefix $(bam_dir)/, $(addprefix lippold_, spy1.bam mez2.bam elsidron2.bam den8.bam den4.bam $(published_bams)))
exome_bams := $(addprefix $(bam_dir)/, $(addprefix exome_, spy1.bam mez2.bam elsidron1.bam den8.bam den4.bam $(published_bams)))

test_bams := $(bam_dir)/control_vindija.bam $(bam_dir)/control_stuttgart.bam

# pileup files
full_pileups := $(addprefix $(pileup_dir)/, $(addprefix full_, spy1.txt.gz mez2.txt.gz den8.txt.gz den4.txt.gz S_French-1.txt.gz))

# VCF files
published_vcfs := $(subst .bam,.vcf.gz, $(published_bams))
full_arch_vcfs    := $(addprefix $(vcf_dir)/, $(addprefix full_, shotgun_spy1.vcf.gz shotgun_mez2.vcf.gz spy1.vcf.gz mez2.vcf.gz den8.vcf.gz den4.vcf.gz))
lippold_arch_vcfs := $(addprefix $(vcf_dir)/, $(addprefix lippold_, spy1.vcf.gz mez2.vcf.gz elsidron2.vcf.gz den8.vcf.gz den4.vcf.gz))
exome_arch_vcfs   := $(addprefix $(vcf_dir)/, $(addprefix exome_, spy1.vcf.gz mez2.vcf.gz elsidron1.vcf.gz den8.vcf.gz den4.vcf.gz))
full_highcov_vcfs     := $(addprefix $(vcf_dir)/, $(addprefix full_, $(published_vcfs)))
lippold_highcov_vcfs  := $(addprefix $(vcf_dir)/, $(addprefix lippold_, $(published_vcfs)))
exome_highcov_vcfs    := $(addprefix $(vcf_dir)/, $(addprefix exome_, $(published_vcfs)))

full_vcf := $(vcf_dir)/full_highcov.vcf.gz
lippold_vcf := $(vcf_dir)/lippold_highcov.vcf.gz
exome_vcf := $(vcf_dir)/exome_highcov.vcf.gz

test_vcfs := $(vcf_dir)/test_gt.vcf.gz

# FASTA files
archaic_fastas := $(addprefix archaic_,full.fa lippold.fa exome.fa)
neand_fastas := $(addprefix neand_,full.fa lippold.fa exome.fa)
den_fastas := $(addprefix den_,full.fa lippold.fa exome.fa)
modern_fastas := $(addprefix modern_,full.fa lippold.fa exome.fa)
present_fastas := $(addprefix present_,full.fa lippold.fa exome.fa)
all_fastas := $(addprefix $(fasta_dir)/all_,$(archaic_fastas) $(neand_fastas) $(den_fastas) $(modern_fastas) $(present_fastas))
tv_fastas := $(addprefix $(fasta_dir)/tv_,$(archaic_fastas) $(neand_fastas) $(den_fastas) $(modern_fastas) $(present_fastas))
fastas := $(all_fastas) $(tv_fastas)

# scripts
bam_sample := /mnt/expressions/mp/bam-sample/bam-sample.py
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
	@echo -e "\tmake clean        -- delete all generated output file"

init: $(dirs) $(full_bed) $(lippold_bed) $(exome_bed) $(full_sites) $(lippold_sites) $(exome_sites)

bam: $(full_bams) $(lippold_bams) $(exome_bams) $(test_bams)

pileup: $(full_pileups)

vcf: $(full_vcf) $(lippold_vcf) $(exome_vcf) $(full_arch_vcfs) $(lippold_arch_vcfs) $(exome_arch_vcfs) $(test_vcfs)

fasta: $(fastas)

diagnostics:
	bams=`cd $(bam_dir); ls *elsidron* *ustishim* *den* *spy* *mez* | grep 'bam$$' | xargs realpath`; \
	mkdir -p $(fig_dir)/damage; \
	cd $(fig_dir)/damage; \
	for b in $$bams; do \
	    qsub -V -cwd -j y -l virtual_free=500M,h_vmem=500M,class=cow /home/mmeyer/perlscripts/solexa/analysis/substitution_patterns.pl $$b; \
	done



#
# BAM processing
#
# $(bam_dir)/%_denisova8sub.bam: $(coord_dir)/capture_%.bed $(bam_dir)/%_denisova8.bam $(bam_dir)/%_comb_neand.bam
# 	cov_neander=$(shell bedtools coverage -a $< -b $(word 3, $^) -d | awk '{sum+=$$5} END { print sum/NR}'); \
# 	cov_denisova=$(shell bedtools coverage -a $< -b $(word 2, $^) -d | awk '{sum+=$$5} END { print sum/NR}'); \
# 	samtools view -b -s `echo $$cov_neander $$cov_denisova | awk '{print $$1/$$2}'` $(word 2, $^) > $@
# 	samtools index $@

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

$(tmp_dir)/den8.bam:
	$(split_and_merge) den8 /mnt/ngs_data/180503_D00829_0138_BCC49NANXX_PEdi_SN_EE_BN_MG/Bustard/BWA/proc1/s_7_sequence_ancient_hg19_evan.bam input/20190207_Ychromosome_Denisova8.txt
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 den8/den8.bam
	mv $(tmp_dir)/den8.uniq.L35MQ25.bam $@
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
	samtools view -h -b /mnt/expressions/mateja/Late_Neandertals/Final_complete_dataset/L35MQ25_per_individual/Mezmaiskaya2_final.L35MQ25.bam Y -o $@
	samtools index $@

$(tmp_dir)/shotgun_spy1.bam:
	samtools view -h -b /mnt/expressions/mateja/Late_Neandertals/Final_complete_dataset/L35MQ25_per_individual/Spy_final.L35MQ25.bam Y -o $@
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

#
# VCF processing
#

$(vcf_dir)/full_highcov.vcf.gz: $(vcf_dir)/full_chimp.vcf.gz $(full_highcov_vcfs)
	bcftools merge $^ | bcftools annotate -x INFO | bcftools view -M 2 -Oz -o $@.all
	bedtools intersect -header -a $@.all -b $(coord_dir)/capture_full.bed | bgzip -c > $@; rm $@.all
	tabix $@

$(vcf_dir)/lippold_highcov.vcf.gz: $(vcf_dir)/lippold_chimp.vcf.gz $(lippold_highcov_vcfs)
	bcftools merge $^ | bcftools annotate -x INFO | bcftools view -M 2 -Oz -o $@.all
	bedtools intersect -header -a $@.all -b $(coord_dir)/capture_lippold.bed | bgzip -c > $@; rm $@.all
	tabix $@

$(vcf_dir)/exome_highcov.vcf.gz: $(vcf_dir)/exome_chimp.vcf.gz $(exome_highcov_vcfs)
	bcftools merge $^ |  bcftools annotate -x INFO | bcftools view -M 2 -Oz -o $@.all
	bedtools intersect -header -a $@.all -b $(coord_dir)/capture_exome.bed | bgzip -c > $@; rm $@.all
	tabix $@

$(vcf_dir)/full_mez2_snpad.vcf.gz:
	cp /mnt/expressions/Janet/YChr/Mez2/chrY.mq25.map50.vcf $(tmp_dir)
	bgzip -c $(tmp_dir)/chrY.mq25.map50.vcf > $@
	tabix $@

# generate genotypes from the Chimp reference genome
$(vcf_dir)/%_chimp.vcf.gz: $(coord_dir)/capture_%.pos
	$(src_dir)/chimp_vcf.sh $< $(basename $@)
	bgzip $(basename $@)
	tabix $@

# genotype samples by consensus calling
$(vcf_dir)/%.vcf.gz: $(bam_dir)/%.bam
	name="$(shell echo $(basename $(notdir $<)) | sed 's/^[a-z]*_//')"; \
	if [[ "$$name" =~ (den|spy|mez|elsidron) ]]; then mincov=1; else mincov=4; fi; \
	$(bam_sample) --bam $< --ref $(ref_genome) --format vcf \
	    --strategy consensus --mincov $$mincov --minbq 20 --minmq 25 \
	    --sample-name $$name --output $(basename $(basename $@))
	bgzip $(basename $@)
	tabix $@

# testing A00 VCF file for comparing bam-sample and bcftools calls
$(vcf_dir)/test_gt.vcf.gz: $(vcf_dir)/full_a00.vcf.gz $(vcf_dir)/full_den8.vcf.gz $(vcf_dir)/full_ustishim.vcf.gz
	samtools mpileup -B -t DP -Q 20 -q 25 -u -f $(ref_genome) $(bam_dir)/full_a00.bam \
		| bcftools call --ploidy 1 -m -V indels -Oz -o $(tmp_dir)/bcftools_a00.vcf.gz
	tabix $(tmp_dir)/bcftools_a00.vcf.gz
	samtools mpileup -B -t DP -Q 20 -q 25 -u -f $(ref_genome) $(bam_dir)/full_den8.bam \
		| bcftools call --ploidy 1 -m -V indels -Oz -o $(tmp_dir)/bcftools_den8.vcf.gz
	tabix $(tmp_dir)/bcftools_den8.vcf.gz
	samtools mpileup -B -t DP -Q 20 -q 30 -u -f $(ref_genome) $(bam_dir)/full_ustishim.bam \
		| bcftools call --ploidy 1 -m -V indels -Oz -o $(tmp_dir)/bcftools_ustishim.vcf.gz
	tabix $(tmp_dir)/bcftools_ustishim.vcf.gz
	bcftools merge $(vcf_dir)/full_a00.vcf.gz $(tmp_dir)/bcftools_a00.vcf.gz $(vcf_dir)/full_den8.vcf.gz $(tmp_dir)/bcftools_den8.vcf.gz $(vcf_dir)/full_ustishim.vcf.gz $(tmp_dir)/bcftools_ustishim.vcf.gz \
	    | bcftools view -v snps \
	    | bcftools annotate -x INFO,FORMAT/PL \
	    | bcftools reheader -s <(echo -e "consensus_a00\nbcftools_a00\nconsensus_den8\nbcftools_den8\nconsensus_ustishim\nbcftools_ustishim\n"| cat) \
	    | bgzip -c > $@
	tabix $@



#
# FASTA alignments for BEAST analyses
#
# archaics := spy1 mez2 neand den8 kk1 mota bichon loschbour ustishim elsidron1 elsidron2
archaics := neand den8 elsidron1 elsidron2
exclude := spy1 mez2 chimp kk1 mota bichon loschbour chimp S_BedouinB-1 S_Punjabi-1 S_Turkish-1 S_Burmese-1 S_Saami-2 S_Thai-1


$(fasta_dir)/all_archaic_%.fa: $(vcf_dir)/merged_%.vcf.gz
	python $(src_dir)/vcf_to_fasta.py --vcf $< --fasta $@ --exclude $(exclude) --variable

$(fasta_dir)/all_neand_%.fa: $(vcf_dir)/merged_%.vcf.gz
	python $(src_dir)/vcf_to_fasta.py --vcf $< --fasta $@ --exclude den8 $(exclude) --variable

$(fasta_dir)/all_den_%.fa: $(vcf_dir)/merged_%.vcf.gz
	python $(src_dir)/vcf_to_fasta.py --vcf $< --fasta $@ --exclude neand $(exclude) --variable

$(fasta_dir)/all_modern_%.fa: $(vcf_dir)/merged_%.vcf.gz
	python $(src_dir)/vcf_to_fasta.py --vcf $< --fasta $@ --exclude $(exclude) $(archaics) --variable

$(fasta_dir)/all_present_%.fa: $(vcf_dir)/merged_%.vcf.gz
	python $(src_dir)/vcf_to_fasta.py --vcf $< --fasta $@ --exclude ustishim $(exclude) $(archaics) --variable


$(fasta_dir)/tv_archaic_%.fa: $(vcf_dir)/merged_%.vcf.gz
	python $(src_dir)/vcf_to_fasta.py --vcf $< --fasta $@ --exclude $(exclude) --variable --tvonly

$(fasta_dir)/tv_neand_%.fa: $(vcf_dir)/merged_%.vcf.gz
	python $(src_dir)/vcf_to_fasta.py --vcf $< --fasta $@ --exclude den8 $(exclude) --variable --tvonly

$(fasta_dir)/tv_den_%.fa: $(vcf_dir)/merged_%.vcf.gz
	python $(src_dir)/vcf_to_fasta.py --vcf $< --fasta $@ --exclude neand $(exclude) --variable --tvonly

$(fasta_dir)/tv_modern_%.fa: $(vcf_dir)/merged_%.vcf.gz
	python $(src_dir)/vcf_to_fasta.py --vcf $< --fasta $@ --exclude $(exclude) $(archaics) --variable --tvonly

$(fasta_dir)/tv_present_%.fa: $(vcf_dir)/merged_%.vcf.gz
	python $(src_dir)/vcf_to_fasta.py --vcf $< --fasta $@ --exclude ustishim $(exclude) $(archaics) --variable --tvonly


#
# generate pileup files for read-based contamination estimates
#
$(pileup_dir)/%.txt.gz: $(bam_dir)/%.bam
	$(bam_sample) --bam $< --ref $(ref_genome) --format pileup \
	    --minbq 20 --minmq 25 \
	    --output $(basename $(basename $@))
	bgzip $(basename $@)



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
