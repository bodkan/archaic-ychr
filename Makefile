SHELL := /bin/bash

doc_dir := doc

bam_dir := bam
tmp_dir := tmp
input_dir := input
figures_dir := figures
output_dir := output
vcf_dir := vcf
data_dirs := $(bam_dir) $(vcf_dir) $(figures_dir) $(input_dir) $(output_dir) $(tmp_dir)

targets_bed := $(input_dir)/target_regions.bed

sidron_bam      := $(bam_dir)/sidron_ontarget.bam
exome_sidron_bam      := $(bam_dir)/exome_sidron_ontarget.bam
den8_bam        := $(bam_dir)/den8_ontarget.bam
deam_den8_bam   := $(bam_dir)/deam_den8_ontarget.bam
a00_bam := $(bam_dir)/a00_ontarget.bam

sidron_vcf := $(vcf_dir)/sidron_ontarget.vcf.gz
a00_vcf := $(vcf_dir)/a00_ontarget.vcf.gz

nb_sidron_processing := $(doc_dir)/processing_of_El_Sidron_data.ipynb
nb_den8_processing := $(doc_dir)/processing_of_Denisova_8_data.ipynb
nb_ancient_features := $(doc_dir)/aDNA_features_analysis.ipynb
nb_coverage_analysis := $(doc_dir)/capture_efficiency_and_coverage__Python.ipynb

ref_genome := /mnt/solexa/Genomes/hg19_evan/whole_genome.fa


.PHONY: default init clean clean_all


default:
	@echo -e "Usage:"
	@echo -e "\tmake init              -- create all necessary directories"
	@echo -e "\tmake bam_processing    -- process all BAM files for the analysis"
	@echo -e "\tmake genotypes         -- run genotyping on all processed BAM files"
	@echo -e "\tmake ancient_features  -- analyze patterns of ancient DNA damage"
	@echo -e "\tmake coverage_analysis -- analyze patterns of ancient DNA damage"


init: $(data_dirs)

process_bams: $(data_dirs) $(sidron_bam) $(den8_bam) $(deam_den8_bam) $(exome_sidron_bam) $(a00_bam)

genotypes: $(data_dirs) $(sidron_vcf) $(a00_vcf)

ancient_features: $(data_dirs)
	jupyter nbconvert $(nb_ancient_features) --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $(nb_ancient_features)

coverage_analysis: $(data_dirs)
	jupyter nbconvert $(nb_coverage_analysis) --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $(nb_coverage_analysis)


$(sidron_bam): $(targets_bed)
	jupyter nbconvert $(nb_sidron_processing) --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $(nb_sidron_processing)

$(den8_bam) $(den8_deam_bam): $(targets_bed)
	jupyter nbconvert $(nb_den8_processing) --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $(nb_den8_processing)

$(exome_sidron_bam):
	cd $(tmp_dir); \
	curl -O http://cdna.eva.mpg.de/neandertal/exomes/BAM/Sidron_exome_hg19_1000g_LowQualDeamination.md.bam; \
	cd ..; \
	bedtools intersect -a $(tmp_dir)/Sidron_exome_hg19_1000g_LowQualDeamination.md.bam -b $(targets_bed) \
		> $@; \
	samtools index $@

$(a00_bam): $(tmp_dir)/A00.bam
	bedtools intersect -a $< -b $(targets_bed) \
		> $@; \
	samtools index $@

$(tmp_dir)/A00.bam:
	cp /mnt/genotyping/sendru/Chiara/validation/A00.bam $@

$(sidron_vcf): $(sidron_bam)
	samtools mpileup -A -Q 20 -u -f $(ref_genome) $< \
		| bcftools call --ploidy 1 -m -V indels -Oz \
		| bcftools reheader -s <(echo -e "ElSidron"| cat) \
		> $@; \
	tabix $@

$(a00_vcf): $(a00_bam)
	samtools mpileup -A -Q 20 -u -f $(ref_genome) $< \
		| bcftools call --ploidy 1 -m -V indels -Oz \
		| bcftools reheader -s <(echo -e "A00"| cat) \
		> $@; \
	tabix $@

$(targets_bed):
	cp /mnt/454/Carbon_beast_QM/QF_chrY_region.bed $@

$(data_dirs):
	mkdir -p $@



clean:
	rm -rf $(vcf_dir) $(figures_dir) $(input_dir) $(output_dir)

clean_all:
	rm -rf $(data_dirs)
