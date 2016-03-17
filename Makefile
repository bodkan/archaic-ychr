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
humans_bams      := $(wildcard /mnt/scratch/basti/HGDP_chrY_data/raw_data_submission/*.bam)

all_bams := $(mez2_bam) $(spy_bam) $(sidron_bam) $(sidron_dq_bam) $(exome_sidron_bam) $(den8_bam) $(deam_den8_bam) $(den4_bam) $(deam_den4_bam) $(a00_1_bam) $(a00_2_bam)
all_bais := $(addsuffix .bai,$(all_bams))

#
# VCF files
#
chimp_vcf      := $(vcf_dir)/chimp_ontarget.vcf.gz

mez2_vcf       := $(vcf_dir)/mez2_ontarget.vcf.gz
spy_vcf        := $(vcf_dir)/spy_ontarget.vcf.gz
sidron_vcf     := $(vcf_dir)/sidron_ontarget.vcf.gz
sidron_dq_vcf  := $(vcf_dir)/sidron_dq_ontarget.vcf.gz
sidron_maj_vcf := $(vcf_dir)/sidron_maj_ontarget.vcf.gz

den8_vcf       := $(vcf_dir)/den8_ontarget.vcf.gz
deam_den8_vcf  := $(vcf_dir)/deam_den8_ontarget.vcf.gz

a00_1_vcf      := $(vcf_dir)/a00_1_ontarget.vcf.gz
a00_2_vcf      := $(vcf_dir)/a00_2_ontarget.vcf.gz
humans_vcf     := $(vcf_dir)/humans_ontarget.vcf.gz

merged_all_vcf := $(vcf_dir)/merged_all_ontarget.vcf.gz
merged_var_vcf := $(vcf_dir)/merged_var_ontarget.vcf.gz

all_vcfs := $(chimp_vcf) $(mez2_vcf) $(spy_vcf) $(sidron_vcf) $(sidron_dq_vcf) $(sidron_maj_vcf) $(den8_vcf) $(deam_den8_vcf) $(a00_1_vcf) $(a00_2_vcf) $(humans_vcf)
all_tbis := $(addsuffix .tbi,$(all_vcfs))

#
# FASTA files
#
chimp_nea_humans_all_sites_fasta    := $(fasta_dir)/all_sites__chimp_nea_humans.fa
chimp_nea_humans_var_sites_fasta    := $(fasta_dir)/var_sites__chimp_nea_humans.fa
chimp_sidron_humans_all_sites_fasta := $(fasta_dir)/all_sites__chimp_sidron_humans.fa
chimp_sidron_humans_var_sites_fasta := $(fasta_dir)/var_sites__chimp_sidron_humans.fa
sidron_humans_all_sites_fasta       := $(fasta_dir)/all_sites__sidron_humans.fa
sidron_humans_var_sites_fasta       := $(fasta_dir)/var_sites__sidron_humans.fa
humans_all_sites_fasta			    := $(fasta_dir)/all_sites__humans.fa

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
bam_sample := ~/devel/bam-utils/bam-sample.py

#
# other files
#
ref_genome := /mnt/solexa/Genomes/hg19_evan/whole_genome.fa

targets_bed := $(input_dir)/target_regions.bed
target_sites := $(input_dir)/target_sites.bed

sample_info = $(input_dir)/sample_info.tsv

# subset of 623 Y chromosome data from Lippold et al.
humans_subset := HGDP00001 HGDP00099 HGDP00449 HGDP00511 HGDP00540 HGDP00608 HGDP00703 HGDP00786

.PHONY: default init clean clean_all


default:
	@echo -e "Usage:"
	@echo -e "\tmake init              -- create all necessary directories"
	@echo -e "\tmake bams              -- process all BAM files for the analysis"
	@echo -e "\tmake genotypes         -- run genotyping on all processed BAM files"
	@echo -e "\tmake ancient_features  -- analyze patterns of ancient DNA damage"
	@echo -e "\tmake coverage_analysis -- analyze patterns of ancient DNA damage"
	@echo -e "\tmake alignments        -- generate FASTA alignments from VCF files"
	@echo -e "\tmake clean             -- delete all generated output file"


init: $(data_dirs)

bams: $(data_dirs) $(all_bams) $(all_bais)

genotypes: $(data_dirs) $(merged_all_vcf) $(merged_var_vcf) $(merged_all_vcf).tbi $(merged_var_vcf).tbi

alignments: $(data_dirs) $(chimp_nea_humans_all_sites_fasta) $(chimp_nea_humans_var_sites_fasta) $(chimp_sidron_humans_all_sites_fasta) $(chimp_sidron_humans_var_sites_fasta) $(sidron_humans_all_sites_fasta) $(sidron_humans_var_sites_fasta) $(humans_all_sites_fasta)

ancient_features: $(data_dirs)
	jupyter nbconvert $(nb_ancient_features) --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $(nb_ancient_features)

coverage_analysis: $(data_dirs)
	jupyter nbconvert $(nb_coverage_analysis) --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $(nb_coverage_analysis)


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

$(exome_sidron_bam):
	cd $(tmp_dir); \
	curl -O http://cdna.eva.mpg.de/neandertal/exomes/BAM/Sidron_exome_hg19_1000g_LowQualDeamination.md.bam; \
	cd ..; \
	bedtools intersect -a $(tmp_dir)/Sidron_exome_hg19_1000g_LowQualDeamination.md.bam -b $(targets_bed) \
		> $@; \
	samtools index $@

$(a00_1_bam): $(tmp_dir)/GRC13292545.chrY.bam
	bedtools intersect -a $< -b $(targets_bed) \
		> $@; \
	samtools index $@

$(a00_2_bam):  $(tmp_dir)/GRC13292546.chrY.bam
	bedtools intersect -a $< -b $(targets_bed) \
		> $@; \
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

$(sidron_maj_vcf): $(sidron_bam)
	python3 $(bam_sample) \
		--bam $(sidron_bam) --bed $(targets_bed) \
		--ref $(ref_genome) --format VCF --sample-name ElSidronMaj \
		--strand-check USER_term5 --method majority \
		--minbq 20 --mincov 3 \
	| bgzip \
	> $@

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

$(humans_vcf): $(humans_bams)
	samtools mpileup -l $(targets_bed) -A -Q 20 -u -f $(ref_genome) $^ \
		|  bcftools call --ploidy 1 -m -V indels -Oz -o $@

$(merged_all_vcf): $(all_vcfs) $(all_tbis)
	bcftools merge -m all $(all_vcfs)  \
		| bcftools view -M2 \
		| bcftools annotate -x INFO,FORMAT/PL -Oz -o $@

$(merged_var_vcf): $(all_vcfs) $(all_tbis)
	bcftools merge -m all $(mez2_vcf) $(spy_vcf) $(sidron_vcf) $(sidron_maj_vcf) $(den8_vcf) $(deam_den8_vcf) $(a00_1_vcf) $(a00_2_vcf) $(humans_vcf) \
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
sample_ids := ElSidron ElSidronMaj A00_1 A00_2 $(humans_subset)

$(chimp_nea_humans_all_sites_fasta): $(merged_all_vcf)
	python $(src_dir)/vcf_to_fasta.py --vcf-file $< --fasta-file $@ --chrom Y --sample-names Chimp Mez2 Spy $(sample_ids)

$(chimp_nea_humans_var_sites_fasta): $(merged_var_vcf)
	python $(src_dir)/vcf_to_fasta.py --vcf-file $< --fasta-file $@ --chrom Y --sample-names Chimp Mez2 Spy $(sample_ids)


$(chimp_sidron_humans_all_sites_fasta): $(merged_all_vcf)
	python $(src_dir)/vcf_to_fasta.py --vcf-file $< --fasta-file $@ --chrom Y --sample-names Chimp $(sample_ids)

$(chimp_sidron_humans_var_sites_fasta): $(merged_var_vcf)
	python $(src_dir)/vcf_to_fasta.py --vcf-file $< --fasta-file $@ --chrom Y --sample-names Chimp $(sample_ids)


$(sidron_humans_all_sites_fasta): $(merged_all_vcf)
	python $(src_dir)/vcf_to_fasta.py --vcf-file $< --fasta-file $@ --chrom Y --sample-names $(sample_ids)

$(sidron_humans_var_sites_fasta): $(merged_var_vcf)
	python $(src_dir)/vcf_to_fasta.py --vcf-file $< --fasta-file $@ --chrom Y --sample-names $(sample_ids)


$(humans_all_sites_fasta): $(merged_all_vcf)
	python $(src_dir)/vcf_to_fasta.py --vcf-file $< --fasta-file $@ --chrom Y --sample-names A00_1 A00_2 $(humans_subset)


#
# other things
#
$(targets_bed):
	cp /mnt/454/Carbon_beast_QM/QF_chrY_region.bed $@

$(target_sites): $(targets_bed)
	python $(src_dir)/sites_in_bed.py --bed-file $< --output-file $@ --format BED

$(sample_info):
	python3 -c "import pandas; df = pandas.read_excel('http://static-content.springer.com/esm/art%3A10.1186%2F2041-2223-5-13/MediaObjects/13323_2014_104_MOESM1_ESM.xlsx', skiprows=6, header=None, parse_cols=[0,1,2]); df.columns = ['name', 'popul', 'region']; df.to_csv('$@', sep='\t', index=False)"



$(data_dirs):
	mkdir -p $@


clean:
	rm -rf $(data_dirs)
