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
sgdp_bams :=  S_French-1.Y.bam S_Sardinian-1.Y.bam S_Han-2.Y.bam S_Dai-2.Y.bam S_Papuan-2.Y.bam S_Karitiana-1.Y.bam S_Dinka-1.Y.bam S_Mbuti-1.Y.bam S_Yoruba-2.Y.bam S_Mandenka-1.Y.bam
published_bams := $(sgdp_bams) ustishim.bam a00_1.bam a00_2.bam
full_bams := $(addprefix $(bam_dir)/, $(addprefix full_, $(published_bams)))
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

damage_patterns:
	bams=$(shell ls $(bam_dir)/*.bam); \
	cd $(fig_dir); \
	for $$bam in $$bams; do \
	    /home/mmeyer/perlscripts/solexa/analysis/substitution_patterns.pl $$bam; \
	done



#
# BAM processing
#
$(bam_dir)/exome_%.bam: $(tmp_dir)/%.bam
	bedtools intersect -a $< -b $(coord_dir)/capture_exome.bed > $@
	samtools index $@

$(bam_dir)/lippold_%.bam: $(tmp_dir)/%.bam
	bedtools intersect -a $< -b $(coord_dir)/capture_lippold.bed > $@
	samtools index $@

$(addprefix $(tmp_dir)/, $(sgdp_bams)):
	cp /mnt/genotyping/sk_pipelines/datasets/Mallick2016_SGDP_Ychromosome/$(notdir $@) $@

$(tmp_dir)/ustishim.bam:
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 /mnt/454/Vindija/high_cov/final_bam/Ust_Ishim/chrY.bam
	mv $(tmp_dir)/chrY.uniq.L35MQ25.bam $@
	samtools index $@

# A00 Y
$(tmp_dir)/a00_1.bam:
	curl -o $@ http://evolbio.ut.ee/chrY/GRC13292545.chrY.bam

$(tmp_dir)/a00_2.bam:
	curl -o $@ http://evolbio.ut.ee/chrY/GRC13292546.chrY.bam

$(tmp_dir)/elsidron1.bam:
	curl http://cdna.eva.mpg.de/neandertal/exomes/BAM/Sidron_exome_hg19_1000g_LowQualDeamination.md.bam -o $@; \
	samtools index $@

$(tmp_dir)/elsidron2.bam:
	$(split_and_merge) elsidron_run1 /mnt/ngs_data/130917_SN7001204_0228_BH06Y0ADXX_R_PEdi_A3207_A3208/Ibis/BWA/s_2-hg19_evan.bam input/elsidron.txt
	$(split_and_merge) elsidron_run2 /mnt/ngs_data/131129_SN7001204_0235_BH72E4ADXX_R_PEdi_A3601_A3605/Bustard/BWA/s_2_sequence_ancient_hg19_evan.bam input/elsidron.txt
	samtools merge $(tmp_dir)/elsidron_run1/*.bam $(tmp_dir)/elsidron_run2/*.bam > $(tmp_dir)/elsidron_both_runs.bam
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 $(tmp_dir)/elsidron_both_runs.bam
	mv $(tmp_dir)/elsidron_both_runs.uniq.L35Q25.bam $@



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
