SHELL := /bin/bash

# directories
data_dir := data
bam_dir := $(data_dir)/bam
fasta_dir := $(data_dir)/fasta
vcf_dir := $(data_dir)/vcf
coord_dir := $(data_dir)/coord
tmp_dir := tmp
dep_dir := dep
fig_dir := fig
src_dir := src
dirs := $(data_dir) $(bam_dir) $(vcf_dir) $(fasta_dir) $(coord_dir) $(fig_dir) $(dep_dir) $(tmp_dir)

# BAM files
sgdp_bams :=  S_French-1.bam S_Sardinian-1.bam S_Han-2.bam S_Dai-2.bam S_Papuan-2.bam S_Karitiana-1.bam S_Dinka-1.bam S_Mbuti-1.bam S_Yoruba-2.bam S_Mandenka-1.bam
published_bams := $(sgdp_bams) ustishim.bam a00_1.bam a00_2.bam kk1.bam mota.bam bichon.bam loschbour.bam
full_bams := $(addprefix $(bam_dir)/, $(addprefix full_, spy1.bam mez2.bam denisova8.bam $(published_bams)))
lippold_bams := $(addprefix $(bam_dir)/, $(addprefix lippold_, elsidron2.bam $(published_bams)))
exome_bams := $(addprefix $(bam_dir)/, $(addprefix exome_, elsidron1.bam $(published_bams)))

# VCF files
all_vcfs := mez2.vcf.gz spy.vcf.gz elsidron.vcf.gz ustishim.vcf.gz a00.vcf.gz
exome_vcfs := $(addprefix $(vcf_dir)/, $(addprefix exome_,$(all_vcfs)))

# FASTA files
fastas := chimp_nea_bteam.fa sidron_bteam.fa bteam.fa
exome_fastas := $(addprefix $(fasta_dir)/exome_,$(fastas))

# scripts
bam_sample := $(dep_dir)/bam-sample/bam-sample
run_nb := $(src_dir)/run_nb.sh
split_and_merge := $(src_dir)/split_and_merge.sh
split_bam := /r1/people/mmeyer/perlscripts/solexa/filework/splitBAM.pl
analyze_bam := /home/mmeyer/perlscripts/solexa/analysis/analyzeBAM.pl

# coordinates
full_bed := $(coord_dir)/capture_full.bed
lippold_bed := $(coord_dir)/capture_lippold.bed
exome_bed := $(coord_dir)/capture_exome.bed
annot_bed := $(coord_dir)/cds.bed $(coord_dir)/phastcons.bed $(coord_dir)/genes.bed $(coord_dir)/pseudogenes.bed
full_sites := $(coord_dir)/sites_full.pos
lippold_sites := $(coord_dir)/sites_lippold.pos
exome_sites := $(coord_dir)/sites_exome.bed

ref_genome := /mnt/solexa/Genomes/hg19_evan/whole_genome.fa

map_filter := /mnt/454/HCNDCAM/Hengs_Alignability_Filter/hs37m_filt35_99.bed.gz


.PHONY: default init bam vcf fasta diagnostics clean



default:
	@echo -e "Usage:"
	@echo -e "\tmake init         -- create all necessary directories"
	@echo -e "\tmake bam          -- process and filter BAM files"
	@echo -e "\tmake vcf          -- run consensus-based genotyping"
	@echo -e "\tmake fasta        -- generate FASTA alignments from VCF files"
	@echo -e "\tmake diagnostics  -- generate diagnostic plots on BAMs"
	@echo -e "\tmake clean        -- delete all generated output file"

init: $(dirs) $(full_bed) $(lippold_bed) $(exome_bed) $(annot_bed) $(bam_sample)

bam: $(dirs) $(full_bams) $(lippold_bams) $(exome_bams)

vcf: $(dirs) $(full_vcfs) $(lippold_vcfs) $(exome_vcfs)

fasta: $(dirs) $(full_fastas) $(lippold_fastas) $(exome_fastas)

diagnostics:
	@bams=`ls $(bam_dir)/*.bam`; \
	mkdir $(fig_dir)/damage; \
	cd $(fig_dir)/damage; \
	for $$bam in $$bams; do \
	    /home/mmeyer/perlscripts/solexa/analysis/substitution_patterns.pl ../../$${bam}; \
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
	cp /mnt/genotyping/sk_pipelines/datasets/Mallick2016_SGDP_Ychromosome/$(basename $(notdir $@)).Y.bam $@

# A00 Y
$(tmp_dir)/a00_1.bam:
	cd $(tmp_dir); curl http://evolbio.ut.ee/chrY/GRC13292545.chrY.bam
	bedtools intersect -a $(tmp_dir)/GRC13292545.chrY.bam -b $(map_filter) > $@

$(tmp_dir)/a00_2.bam:
	cd $(tmp_dir); curl http://evolbio.ut.ee/chrY/GRC13292546.chrY.bam
	bedtools intersect -a $(tmp_dir)/GRC13292546.chrY.bam -b $(map_filter) > $@

$(tmp_dir)/elsidron1.bam:
	curl http://cdna.eva.mpg.de/neandertal/exomes/BAM/Sidron_exome_hg19_1000g_LowQualDeamination.md.bam -o $@; \
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 $(notdir $@)
	bedtools intersect -a $(tmp_dir)/elsidron1.uniq.L35MQ25.bam -b $(map_filter) > $@
	samtools index $@

$(tmp_dir)/elsidron2.bam: $(tmp_dir)/elsidron_run1/elsidron_run1.bam $(tmp_dir)/elsidron_run2/elsidron_run2.bam
	samtools merge $(tmp_dir)/elsidron_both_runs.bam $(tmp_dir)/elsidron_run1/elsidron_run1.bam $(tmp_dir)/elsidron_run2/elsidron_run2.bam
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 elsidron_both_runs.bam
	bedtools intersect -a $(tmp_dir)/elsidron_both_runs.uniq.L35MQ25.bam -b $(map_filter) > $@
	samtools index $@

$(tmp_dir)/elsidron_run1/elsidron_run1.bam:
	$(split_and_merge) elsidron_run1 /mnt/ngs_data/130917_SN7001204_0228_BH06Y0ADXX_R_PEdi_A3207_A3208/Ibis/BWA/s_2-hg19_evan.bam input/A2970_A3206_A3208.txt
$(tmp_dir)/elsidron_run2/elsidron_run2.bam:
	$(split_and_merge) elsidron_run2 /mnt/ngs_data/131129_SN7001204_0235_BH72E4ADXX_R_PEdi_A3601_A3605/Bustard/BWA/s_2_sequence_ancient_hg19_evan.bam input/A2970_A3206_A3208.txt

$(tmp_dir)/ustishim.bam:
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 /mnt/454/Vindija/high_cov/final_bam/Ust_Ishim/chrY.bam
	bedtools intersect -a $(tmp_dir)/chrY.uniq.L35MQ25.bam -b $(map_filter) > $@
	samtools index $@

$(tmp_dir)/kk1.bam:
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 /mnt/expressions/mp/Archive/y-selection/tmp/KK1_sort_rmdup_merge_IR_q30_mapDamage.bam
	bedtools intersect -a $(tmp_dir)/KK1_sort_rmdup_merge_IR_q30_mapDamage.uniq.L35MQ25.bam -b $(map_filter) > $@
	samtools index $@

$(tmp_dir)/bichon.bam:
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 /mnt/expressions/mp/Archive/y-selection/tmp/Bichon.sort.rmdup.IR.q30.mapDamage.bam
	bedtools intersect -a $(tmp_dir)/Bichon.sort.rmdup.IR.q30.mapDamage.uniq.L35MQ25.bam -b $(map_filter) > $@
	samtools index $@

$(tmp_dir)/mota.bam:
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 /mnt/expressions/mp/Archive/y-selection/tmp/GB20_sort_merge_dedup_l30_IR_q30_mapDamage.bam
	bedtools intersect -a $(tmp_dir)/GB20_sort_merge_dedup_l30_IR_q30_mapDamage.uniq.L35MQ25.bam -b $(map_filter) > $@
	samtools index $@

$(tmp_dir)/loschbour.bam:
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 /mnt/expressions/mp/Archive/y-selection/tmp/Loschbour.hg19_1000g.bam
	bedtools intersect -a $(tmp_dir)/Loschbour.hg19_1000g.uniq.L35MQ25.bam -b $(map_filter) > $@
	samtools index $@

$(tmp_dir)/denisova8.bam:
	$(split_and_merge) denisova8 /mnt/ngs_data/180503_D00829_0138_BCC49NANXX_PEdi_SN_EE_BN_MG/Bustard/BWA/proc1/s_7_sequence_ancient_hg19_evan.bam input/20190207_Ychromosome_Denisova8.txt
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 denisova8/denisova8.bam
	bedtools intersect -a $(tmp_dir)/denisova8.uniq.L35MQ25.bam -b $(map_filter) > $@
	samtools index $@

$(tmp_dir)/spy1.bam:
	$(split_and_merge) spy1 /mnt/ngs_data/170825_D00829_0064_AHNVL5BCXY_R_PEdi_F5281_F5282/Bustard/BWA/proc1/s_2_sequence_ancient_hg19_evan.bam input/20190207_Ychromosome_Spy1.txt
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 spy1/spy1.bam
	bedtools intersect -a $(tmp_dir)/spy1.uniq.L35MQ25.bam -b $(map_filter) > $@
	samtools index $@

$(tmp_dir)/mez2.bam:
	$(split_and_merge) mez2 /mnt/ngs_data/170825_D00829_0064_AHNVL5BCXY_R_PEdi_F5281_F5282/Bustard/BWA/proc1/s_2_sequence_ancient_hg19_evan.bam input/20190207_Ychromosome_Mez2.txt
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 mez2/mez2.bam
	bedtools intersect -a $(tmp_dir)/mez2.uniq.L35MQ25.bam -b $(map_filter) > $@
	samtools index $@



#
# VCF processing
#



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
$(lippold_bed):
	cp input/basti_design.bed $@

# Y chromosome capture regions designed by Qiaomei (~6Mb)
$(full_bed):
	perl -lane 'print $$F[0] . "\t" . $$F[1] . "\t" . $$F[2]' input/Y.filt35_50_SRepeat_100.bed > $@

# Y chromosome exome capture regions
$(exome_bed):
	wget http://www.cell.com/cms/attachment/2052899616/2060015784/mmc2.zip
	unzip mmc2.zip; rm mmc2.zip
	mv ajhg2064mmc2_V1.txt $@

# functional annotation coordinates
$(annot_bed):
	$(run_nb) notebooks/annotated_regions.ipynb

# sites within BED files
$(coord_dir)/sites_%.pos: $(coord_dir)/capture_%.bed
	$(src_dir)/sites_in_bed.py --bed $< --format pos --output $@


#
# other dependencies
#
$(bam_sample):
	git clone https://github.com/bodkan/bam-sample $(dep_dir)/bam-sample

$(dirs):
	mkdir -p $@

clean:
	rm -rf $(dirs)
