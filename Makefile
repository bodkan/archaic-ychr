SHELL := /bin/bash

#
# directories
#
doc_dir := doc
bam_dir := bam
tmp_dir := tmp
input_dir := input
figures_dir := figures
fasta_dir := fasta
vcf_dir := vcf
src_dir := src
data_dirs := $(bam_dir) $(vcf_dir) $(figures_dir) $(input_dir) $(fasta_dir) $(tmp_dir)

# table with information about males from the B-team
bteam_info := $(input_dir)/bteam_info.tsv

# sample IDs of B-team males (used in VCF file names)
bteam_ids := $(shell cut -f1 $(bteam_info) | tr '\n' ' ')
# sample populations of B-team males (used in VCF headers)
bteam_names := $(shell cut -f3 $(bteam_info) | tr '\n' ' ')

#
# BAM files
#
mez2_bam         := $(bam_dir)/mez2_ontarget.bam
spy_bam          := $(bam_dir)/spy_ontarget.bam
sidron_bam       := $(bam_dir)/sidron_ontarget.bam
sidron_dq_bam    := $(bam_dir)/sidron_dq_ontarget.bam
exome_sidron_bam := $(bam_dir)/exome_sidron_ontarget.bam

den8_bam         := $(bam_dir)/den8_ontarget.bam
deam_den8_bam    := $(bam_dir)/deam_den8_ontarget.bam
den4_bam         := $(bam_dir)/den4_ontarget.bam
deam_den4_bam    := $(bam_dir)/deam_den4_ontarget.bam

a00_1_bam        := $(bam_dir)/a00_1_ontarget.bam
a00_2_bam        := $(bam_dir)/a00_2_ontarget.bam

all_bams := $(mez2_bam) $(spy_bam) $(sidron_bam) $(sidron_dq_bam) $(exome_sidron_bam) $(den8_bam) $(deam_den8_bam) $(den4_bam) $(deam_den4_bam) $(a00_1_bam) $(a00_2_bam)
all_bais := $(addsuffix .bai,$(all_bams))

whole_exome_bam := $(tmp_dir)/whole_exome.bam

#
# VCF files
#
chimp_vcf      := $(vcf_dir)/chimp_ontarget.vcf.gz

mez2_vcf       := $(vcf_dir)/mez2_ontarget.vcf.gz
spy_vcf        := $(vcf_dir)/spy_ontarget.vcf.gz
sidron_vcf     := $(vcf_dir)/sidron_ontarget.vcf.gz
sidron_dq_vcf  := $(vcf_dir)/sidron_dq_ontarget.vcf.gz
sidron_q_vcf   := $(vcf_dir)/sidron_q_ontarget.vcf.gz
sidron_maj_vcf := $(vcf_dir)/sidron_maj_ontarget.vcf.gz
exome_sidron_vcf := $(vcf_dir)/exome_sidron_ontarget.vcf.gz

den8_vcf       := $(vcf_dir)/den8_ontarget.vcf.gz
deam_den8_vcf  := $(vcf_dir)/deam_den8_ontarget.vcf.gz

a00_1_vcf      := $(vcf_dir)/a00_1_ontarget.vcf.gz
a00_2_vcf      := $(vcf_dir)/a00_2_ontarget.vcf.gz
bteam_vcfs     := $(addsuffix .vcf.gz,$(addprefix vcf/bteam_,$(bteam_ids)))

merged_all_vcf := $(vcf_dir)/merged_all_ontarget.vcf.gz
merged_var_vcf := $(vcf_dir)/merged_var_ontarget.vcf.gz

all_vcfs := $(chimp_vcf) $(mez2_vcf) $(spy_vcf) $(sidron_vcf) $(sidron_dq_vcf) $(sidron_q_vcf) $(sidron_maj_vcf) $(den8_vcf) $(deam_den8_vcf) $(a00_1_vcf) $(a00_2_vcf) $(bteam_vcfs)
all_tbis := $(addsuffix .tbi,$(all_vcfs))

#
# FASTA files
#
chimp_nea_bteam_all_sites_fasta    := $(fasta_dir)/all_sites__chimp_nea_bteam.fa
chimp_nea_bteam_var_sites_fasta    := $(fasta_dir)/var_sites__chimp_nea_bteam.fa
chimp_sidron_bteam_all_sites_fasta := $(fasta_dir)/all_sites__chimp_sidron_bteam.fa
chimp_sidron_bteam_var_sites_fasta := $(fasta_dir)/var_sites__chimp_sidron_bteam.fa
sidron_bteam_all_sites_fasta       := $(fasta_dir)/all_sites__sidron_bteam.fa
sidron_bteam_var_sites_fasta       := $(fasta_dir)/var_sites__sidron_bteam.fa
bteam_all_sites_fasta			   := $(fasta_dir)/all_sites__bteam.fa

#
# Jupyter notebooks used for processing and analysis
#
nb_sidron_processing    := $(doc_dir)/processing_of_El_Sidron_data.ipynb
nb_den_processing       := $(doc_dir)/processing_of_Denisova_shotgun_data.ipynb
nb_ancient_features     := $(doc_dir)/aDNA_features_analysis.ipynb
nb_coverage_analysis    := $(doc_dir)/capture_efficiency_and_coverage__Python.ipynb
nb_chimpanzee_genotypes := $(doc_dir)/get_chimpanzee_genotypes.ipynb

#
# scripts and binaries
#
bam_sample := ~/devel/bam-utils/bam_sample.py
bam_plotdamage := ~/devel/bam-utils/bam_plotdamage.py
process_bteam_vcf := $(src_dir)/process_bteam_vcf.sh

#
# other files
#
ref_genome := /mnt/solexa/Genomes/hg19_evan/whole_genome.fa

targets_bed := $(input_dir)/target_regions_map35-99.bed
target_sites := $(input_dir)/target_sites.bed
exome_targets_bed := $(input_dir)/exome_target_regions.bed


.PHONY: default init clean clean_all


default:
	@echo -e "Usage:"
	@echo -e "\tmake init              -- create all necessary directories"
	@echo -e "\tmake bams              -- process all BAM files for the analysis"
	@echo -e "\tmake genotypes         -- run genotyping on all processed BAM files"
	@echo -e "\tmake ancient_features  -- analyze patterns of ancient DNA damage"
	@echo -e "\tmake coverage_analysis -- analyze patterns of ancient DNA damage"
	@echo -e "\tmake alignments        -- generate FASTA alignments from VCF files"
	@echo -e "\tmake damage_patterns   -- generate plots with damage patterns"
	@echo -e "\tmake clean             -- delete all generated output file"


init: $(data_dirs)

bams: $(data_dirs) $(all_bams) $(all_bais)

genotypes: $(data_dirs) $(merged_all_vcf) $(merged_var_vcf) $(merged_all_vcf).tbi $(merged_var_vcf).tbi $(exome_sidron_vcf) $(exome_sidron_vcf).tbi

alignments: $(data_dirs) $(sidron_bteam_all_sites_fasta) $(bteam_all_sites_fasta)

ancient_features: $(data_dirs)
	jupyter nbconvert $(nb_ancient_features) --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $(nb_ancient_features)

coverage_analysis: $(data_dirs)
	jupyter nbconvert $(nb_coverage_analysis) --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $(nb_coverage_analysis)

damage_patterns: bam/Sidron.hg19_evan.Y.dq_ontarget.bam $(a00_1_bam) $(sidron_bam) $(exome_sidron_bam)
	@cd $(figures_dir); \
	for bam in $^; do \
		echo $$bam; \
		python3 $(bam_plotdamage) --bam ../$$bam; \
	done

#
# BAM processing
#
$(mez2_bam): $(targets_bed)
	bedtools intersect -a /mnt/expressions/mateja/Late_Neandertals/Final_complete_dataset/Merged_per_individual_L35MQ0/Mezmaiskaya2_final.bam -b $< -sorted \
		> $@

$(spy_bam): $(targets_bed)
	bedtools intersect -a /mnt/expressions/mateja/Late_Neandertals/Final_complete_dataset/Merged_per_individual_L35MQ0/Spy_final.bam -b $< -sorted \
		> $@

$(sidron_bam): $(targets_bed)
	jupyter nbconvert $(nb_sidron_processing) --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $(nb_sidron_processing);

$(sidron_dq_bam): $(sidron_bam)
	/r1/people/gabriel_renaud/scripts/libbam/decrQualDeaminatedDoubleStranded -n 5 $< $@


$(den8_bam) $(deam_den8_bam) $(den4_bam) $(deam_den4_bam): $(nb_den_processing) $(targets_bed)
	jupyter nbconvert $< --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $<; \
	touch $(den8_bam) $(deam_den8_bam) $(den4_bam) $(deam_den4_bam)

$(exome_sidron_bam): $(exome_targets_bed) $(tmp_dir)/whole_exome.bam
	bedtools intersect -a $(tmp_dir)/whole_exome.bam -b $< -sorted \
		> $@

$(a00_1_bam): $(tmp_dir)/GRC13292545.chrY.bam
	bedtools intersect -a $< -b $(targets_bed) \
		> $@; \
	samtools index $@

$(a00_2_bam):  $(tmp_dir)/GRC13292546.chrY.bam
	bedtools intersect -a $< -b $(targets_bed) \
		> $@; \
	samtools index $@

$(tmp_dir)/whole_exome.bam:
	curl http://cdna.eva.mpg.de/neandertal/exomes/BAM/Sidron_exome_hg19_1000g_LowQualDeamination.md.bam -o $@; \
	samtools index $@

$(tmp_dir)/GRC13292545.chrY.bam:
	cd $(tmp_dir); curl -O http://evolbio.ut.ee/chrY/GRC13292545.chrY.bam

$(tmp_dir)/GRC13292546.chrY.bam:
	cd $(tmp_dir); curl -O http://evolbio.ut.ee/chrY/GRC13292546.chrY.bam


#
# VCF processing
#
$(chimp_vcf): $(target_sites) $(nb_chimpanzee_genotypes)
	jupyter nbconvert $(nb_chimpanzee_genotypes) --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $(nb_chimpanzee_genotypes); \
	touch $@

$(sidron_vcf): $(sidron_bam)
	samtools mpileup -l $(targets_bed) -A -Q 20 -u -f $(ref_genome) $< \
		| bcftools call --ploidy 1 -m -V indels -Oz \
		| bcftools reheader -s <(echo -e "ElSidron"| cat) -o $@

$(sidron_dq_vcf): $(sidron_dq_bam)
	samtools mpileup -l $(targets_bed) -A -Q 20 -u -f $(ref_genome) $< \
		| bcftools call --ploidy 1 -m -V indels -Oz \
		| bcftools reheader -s <(echo -e "ElSidronDQ"| cat) -o $@

$(sidron_q_vcf): $(targets_bed)
	 bcftools view -V indels /mnt/454/Carbon_beast_QM/Y_Sidron_TY/1_Extended_VCF/Sidron.hg19_evan.Y.mod.vcf.gz -R $< \
		| grep -v "LowQual" \
		| sed 's/0\/0/0/; s/1\/1/1/; s/\.\/\./\./' \
		| grep -v "0\/1" \
		| bcftools annotate -x INFO,FORMAT -Oz -o $@

$(sidron_maj_vcf): $(sidron_bam)
	python3 $(bam_sample) \
		--bam $(sidron_bam) --bed $(targets_bed) \
		--ref $(ref_genome) --format VCF --sample-name ElSidronMaj \
		--strand-check USER_term5 --method majority \
		--minbq 20 --mincov 3 \
	| bgzip \
	> $@

$(exome_sidron_vcf): $(exome_sidron_bam)
	samtools mpileup -l $(exome_targets_bed) -A -Q 20 -u -f $(ref_genome) $< \
		| bcftools call --ploidy 1 -m -V indels -Oz \
		| bcftools reheader -s <(echo -e "ElSidronExome"| cat) -o $@; \
	tabix $@

$(mez2_vcf): $(target_regions) $(mez2_bam)
	python3 $(bam_sample) \
		--bam $(mez2_bam) --bed $(targets_bed) \
		--ref $(ref_genome) --format VCF --sample-name Mez2 \
		--strand-check non-USER_all --method majority \
	| bgzip \
	> $@

$(spy_vcf): $(target_regions) $(spy_bam)
	python3 $(bam_sample) \
		--bam $(spy_bam) --bed $(targets_bed) \
		--ref $(ref_genome) --format VCF --sample-name Spy \
		--strand-check non-USER_all --method majority \
	| bgzip \
	> $@

$(den8_vcf): $(target_regions) $(den8_bam)
	python3 $(bam_sample) \
		--bam $(den8_bam) --bed $(targets_bed) \
		--ref $(ref_genome) --format VCF --sample-name Den8 \
		--strand-check USER --method majority \
	| bgzip \
	> $@

$(deam_den8_vcf): $(target_regions) $(deam_den8_bam)
	python3 $(bam_sample) \
		--bam $(deam_den8_bam) --bed $(targets_bed) \
		--ref $(ref_genome) --format VCF --sample-name Den8_deam \
		--strand-check USER --method majority \
	| bgzip \
	> $@

$(a00_1_vcf): $(a00_1_bam)
	samtools mpileup -l $(targets_bed) -A -Q 20 -u -f $(ref_genome) $< \
		| bcftools call --ploidy 1 -m -V indels -Oz \
		| bcftools reheader -s <(echo -e "A00_1"| cat) -o $@

$(a00_2_vcf): $(a00_2_bam)
	samtools mpileup -l $(targets_bed) -A -Q 20 -u -f $(ref_genome) $< \
		| bcftools call --ploidy 1 -m -V indels -Oz \
		| bcftools reheader -s <(echo -e "A00_2"| cat) -o $@

$(vcf_dir)/bteam_%.vcf.gz: $(bteam_info) $(targets_bed)
	id=$(subst .vcf.gz,,$(subst bteam_,,$(notdir $@))); \
	$(process_bteam_vcf) $$id $^

$(merged_all_vcf): $(all_vcfs) $(all_tbis)
	bcftools merge -m all $(all_vcfs)  \
		| bcftools view -M2 \
		| bcftools annotate -x INFO,FORMAT/PL -Oz -o $@

$(merged_var_vcf): $(all_vcfs) $(all_tbis)
	bcftools merge -m all $(mez2_vcf) $(spy_vcf) $(sidron_vcf) $(sidron_dq_vcf) $(sidron_maj_vcf) $(sidron_q_vcf) $(den8_vcf) $(deam_den8_vcf) $(a00_1_vcf) $(a00_2_vcf) $(bteam_vcfs) \
		| bcftools view -m2 -M2 \
		| bcftools annotate -x INFO,FORMAT/PL -Oz -o $@_tmp; \
	bcftools view $(chimp_vcf) -R $@_tmp -Oz -o $(chimp_vcf)_subset; \
	tabix $@_tmp; tabix $(chimp_vcf)_subset; \
	bcftools merge -m all $(chimp_vcf)_subset $@_tmp \
		| bcftools view -m2 -M2 -Oz -o $@; \
	rm $@_tmp $(chimp_vcf)_subset $@_tmp.tbi $(chimp_vcf)_subset.tbi

$(bam_dir)/%.bai: $(bam_dir)/%
	samtools index $<

$(vcf_dir)/%.vcf.gz.tbi: $(vcf_dir)/%.vcf.gz
	tabix -f $<


#
# FASTA generation
#
sample_ids := ElSidronDQ A00_1 A00_2 $(bteam_names)

$(chimp_nea_bteam_all_sites_fasta): $(merged_all_vcf)
	python $(src_dir)/vcf_to_fasta.py --vcf-file $< --fasta-file $@ --chrom Y --sample-names Chimp Mez2 Spy $(sample_ids)

$(chimp_nea_bteam_var_sites_fasta): $(merged_var_vcf)
	python $(src_dir)/vcf_to_fasta.py --vcf-file $< --fasta-file $@ --chrom Y --sample-names Chimp Mez2 Spy $(sample_ids)


$(chimp_sidron_bteam_all_sites_fasta): $(merged_all_vcf)
	python $(src_dir)/vcf_to_fasta.py --vcf-file $< --fasta-file $@ --chrom Y --sample-names Chimp $(sample_ids)

$(chimp_sidron_bteam_var_sites_fasta): $(merged_var_vcf)
	python $(src_dir)/vcf_to_fasta.py --vcf-file $< --fasta-file $@ --chrom Y --sample-names Chimp $(sample_ids)


$(sidron_bteam_all_sites_fasta): $(merged_all_vcf)
	python $(src_dir)/vcf_to_fasta.py --vcf-file $< --fasta-file $@ --chrom Y --sample-names $(sample_ids)

$(sidron_bteam_var_sites_fasta): $(merged_var_vcf)
	python $(src_dir)/vcf_to_fasta.py --vcf-file $< --fasta-file $@ --chrom Y --sample-names $(sample_ids)


$(bteam_all_sites_fasta): $(merged_all_vcf)
	python $(src_dir)/vcf_to_fasta.py --vcf-file $< --fasta-file $@ --chrom Y --sample-names A00_1 A00_2 $(bteam_names)


#
# other things
#
$(targets_bed):
	bedtools intersect -a /mnt/454/Carbon_beast_QM/QF_chrY_region.bed \
		-b /mnt/454/HCNDCAM/Hengs_Alignability_Filter/hs37m_filt35_99.bed.gz \
		| sort -k1,1n -k2,2n \
		> $@

$(exome_targets_bed):
	cp /mnt/scratch/sergi/exome/coordinates/chr/Y_primary_target+tile_margins.MPI.hg19_1000g.bed $@

$(target_sites): $(targets_bed)
	python $(src_dir)/sites_in_bed.py --bed-file $< --output-file $@ --format BED

$(bteam_info):
	grep '^SS.*M$$' /mnt/454/HighCovNeandertalGenome/1_Extended_VCF/Individuals.txt \
		| tr '\t' ' ' \
		| tr -s ' ' \
		| tr ' ' '\t' \
		| cut -f1,2,3 \
		> $@


$(data_dirs):
	mkdir -p $@


clean:
	rm -rf $(data_dirs)
