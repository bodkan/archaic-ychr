SHELL := /bin/bash

doc_dir := doc

bam_dir := bam
tmp_dir := tmp
input_dir := input
figures_dir := figures
output_dir := output
vcf_dir := vcf
src_dir := src
data_dirs := $(bam_dir) $(vcf_dir) $(figures_dir) $(input_dir) $(output_dir) $(tmp_dir)

targets_bed := $(input_dir)/target_regions.bed
target_sites := $(input_dir)/target_sites.pos

sidron_bam      := $(bam_dir)/sidron_ontarget.bam
exome_sidron_bam      := $(bam_dir)/exome_sidron_ontarget.bam
den8_bam        := $(bam_dir)/den8_ontarget.bam
deam_den8_bam   := $(bam_dir)/deam_den8_ontarget.bam
den4_bam        := $(bam_dir)/den4_ontarget.bam
deam_den4_bam   := $(bam_dir)/deam_den4_ontarget.bam
a00_bam := $(bam_dir)/a00_ontarget.bam
hum_623_bams := $(wildcard /mnt/scratch/basti/HGDP_chrY_data/raw_data_submission/*.bam)

all_bams := $(sidron_bam) $(exome_sidron_bam) $(den8_bam) $(deam_den8_bam) $(den4_bam) $(deam_den4_bam) $(a00_bam) $(hum_623_bams)

chimp_vcf := $(vcf_dir)/chimp_ontarget.vcf.gz
sidron_vcf := $(vcf_dir)/sidron_ontarget.vcf.gz
den8_vcf := $(vcf_dir)/den8_ontarget.vcf.gz
a00_vcf := $(vcf_dir)/a00_ontarget.vcf.gz
hum_623_vcf := $(vcf_dir)/hum_623_ontarget.vcf.gz
merged_vcf := $(vcf_dir)/merged_ontarget.vcf.gz

chimp_tbi := $(vcf_dir)/chimp_ontarget.vcf.gz.tbi
sidron_tbi := $(vcf_dir)/sidron_ontarget.vcf.gz.tbi
den8_tbi := $(vcf_dir)/den8_ontarget.vcf.gz.tbi
a00_tbi := $(vcf_dir)/a00_ontarget.vcf.gz.tbi
hum_623_tbi := $(vcf_dir)/hum_623_ontarget.vcf.gz.tbi
merged_tbi := $(vcf_dir)/merged_ontarget.vcf.gz.tbi

all_vcfs := $(chimp_vcf) $(sidron_vcf) $(a00_vcf) $(hum_623_vcf)
all_tbis :=  $(chimp_tbi) $(sidron_tbi) $(a00_tbi) $(hum_623_tbi)

targets_fasta := $(output_dir)/ontarget.fa

nb_sidron_processing := $(doc_dir)/processing_of_El_Sidron_data.ipynb
nb_den_processing := $(doc_dir)/processing_of_Denisova_shotgun_data.ipynb
nb_ancient_features := $(doc_dir)/aDNA_features_analysis.ipynb
nb_coverage_analysis := $(doc_dir)/capture_efficiency_and_coverage__Python.ipynb
nb_chimpanzee_genotypes := $(doc_dir)/get_chimpanzee_genotypes.ipynb

ref_genome := /mnt/solexa/Genomes/hg19_evan/whole_genome.fa


.PHONY: default init clean clean_all


default:
	@echo -e "Usage:"
	@echo -e "\tmake init              -- create all necessary directories"
	@echo -e "\tmake bams              -- process all BAM files for the analysis"
	@echo -e "\tmake genotypes         -- run genotyping on all processed BAM files"
	@echo -e "\tmake ancient_features  -- analyze patterns of ancient DNA damage"
	@echo -e "\tmake coverage_analysis -- analyze patterns of ancient DNA damage"


init: $(data_dirs)

bams: $(data_dirs) $(all_bams)

genotypes: $(data_dirs) $(merged_vcf) $(merged_tbi)

fasta: $(targets_fasta)

ancient_features: $(data_dirs)
	jupyter nbconvert $(nb_ancient_features) --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $(nb_ancient_features)

coverage_analysis: $(data_dirs)
	jupyter nbconvert $(nb_coverage_analysis) --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $(nb_coverage_analysis)



$(sidron_bam): $(targets_bed)
	jupyter nbconvert $(nb_sidron_processing) --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $(nb_sidron_processing)

$(den8_bam) $(deam_den8_bam) $(den4_bam) $(deam_den4_bam): $(nb_den_processing) $(targets_bed)
	jupyter nbconvert $< --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $<

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

$(chimp_vcf): $(targets_sites) $(nb_chimpanzee_genotypes)
	jupyter nbconvert $(nb_chimpanzee_genotypes) --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $(nb_chimpanzee_genotypes); \
	touch $@

$(sidron_vcf): $(sidron_bam)
	samtools mpileup -l $(targets_bed) -A -Q 20 -u -f $(ref_genome) $< \
		| bcftools call --ploidy 1 -m -V indels -Oz \
		| bcftools reheader -s <(echo -e "ElSidron"| cat) -o $@

$(den8_vcf): $(den8_bam)
	samtools mpileup -l $(targets_bed) -A -Q 20 -u -f $(ref_genome) $< \
		| bcftools call --ploidy 1 -m -V indels -Oz \
		| bcftools reheader -s <(echo -e "Denisova8"| cat) -o $@

$(a00_vcf): $(a00_bam)
	samtools mpileup -l $(targets_bed) -A -Q 20 -u -f $(ref_genome) $< \
		| bcftools call --ploidy 1 -m -V indels -Oz \
		| bcftools reheader -s <(echo -e "A00"| cat) -o $@

$(hum_623_vcf): $(hum_623_bams)
	samtools mpileup -l $(targets_bed) -A -Q 20 -u -f $(ref_genome) $^ \
		|  bcftools call --ploidy 1 -m -V indels -Oz -o $@

$(vcf_dir)/%.vcf.gz.tbi: $(vcf_dir)/%.vcf.gz
	tabix -f $<

$(merged_vcf): $(all_vcfs) $(all_tbis)
	bcftools merge -m all $(all_vcfs) \
		| bcftools view -m2 -M2 -Oz -o $@

$(targets_bed):
	cp /mnt/454/Carbon_beast_QM/QF_chrY_region.bed $@

$(target_sites): $(targets_bed)
	python $(src_dir)/sites_in_bed.py --bed-file $< --output-file $@ --format POS

$(targets_fasta): $(merged_vcf)
	python $(src_dir)/vcf_to_fasta.py --vcf-file $< --fasta-file $@ --chrom Y

$(data_dirs):
	mkdir -p $@



clean:
	rm -rf $(vcf_dir) $(figures_dir) $(input_dir) $(output_dir)

clean_all:
	rm -rf $(data_dirs)
