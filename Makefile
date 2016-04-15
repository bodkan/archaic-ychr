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
all_bams := mez2.bam spy.bam sidron.bam a00_1.bam a00_2.bam

lippold_bams := $(addprefix $(bam_dir)/, $(addprefix lippold_,$(all_bams)))
exome_bams   := $(addprefix $(bam_dir)/, $(addprefix exome_,$(all_bams)))
lippold_bais := $(addsuffix .bai,$(lippold_bams))
exome_bais   := $(addsuffix .bai,$(exome_bams))

# Denisova 4/8 BAMS -- not used for actual analysis, but at least to
# know the coverage, number of available sites etc...
denisova_bams := $(addprefix $(bam_dir)/,den8_ontarget.bam deam_den8_ontarget.bam den4_ontarget.bam deam_den4_ontarget.bam)

#
# individual VCF files
#
all_vcfs := chimp.vcf.gz mez2.vcf.gz spy.vcf.gz sidron.vcf.gz a00_1.vcf.gz a00_2.vcf.gz $(addsuffix .vcf.gz,$(addprefix bteam_,$(bteam_ids)))

lippold_vcfs := $(addprefix $(vcf_dir)/, $(addprefix lippold_,$(all_vcfs)))
exome_vcfs   := $(addprefix $(vcf_dir)/, $(addprefix exome_,$(all_vcfs)))
lippold_tbis := $(addsuffix .tbi,$(lippold_vcfs))
exome_tbis   := $(addsuffix .tbi,$(exome_vcfs))

lippold_merged_vcf := $(vcf_dir)/merged_lippold.vcf.gz
exome_merged_vcf   := $(vcf_dir)/merged_exome.vcf.gz
comb_merged_vcf    := $(vcf_dir)/merged_combined.vcf.gz

dpfilt_lippold_merged_vcf := $(vcf_dir)/dpfilt_merged_lippold.vcf.gz
dpfilt_exome_merged_vcf   := $(vcf_dir)/dpfilt_merged_exome.vcf.gz
dpfilt_comb_merged_vcf    := $(vcf_dir)/dpfilt_merged_combined.vcf.gz

#
# FASTA files
#
fastas := chimp_nea_bteam.fa sidron_bteam.fa bteam.fa
lippold_fastas  := $(addprefix $(fasta_dir)/lippold_,$(fastas))
exome_fastas    := $(addprefix $(fasta_dir)/exome_,$(fastas))
combined_fastas := $(addprefix $(fasta_dir)/combined_,$(fastas))

#
# Jupyter notebooks used for processing and analysis
#
nb_sidron_processing    := $(doc_dir)/processing_of_El_Sidron_data.ipynb
nb_ancient_features     := $(doc_dir)/aDNA_features_analysis.ipynb
nb_coverage_analysis    := $(doc_dir)/capture_efficiency_and_coverage__Python.ipynb
nb_chimpanzee_genotypes := $(doc_dir)/get_chimpanzee_genotypes.ipynb
nb_denisova_processing  := $(doc_dir)/processing_of_Denisova_shotgun_data.ipynb

#
# scripts and binaries
#
bam_sample := ~/devel/bam-utils/bam_sample.py
bam_plotdamage := ~/devel/bam-utils/bam_plotdamage.py
process_bteam_vcf := $(src_dir)/process_bteam_vcf.sh
decrease_bquals := /r1/people/gabriel_renaud/scripts/libbam/decrQualDeaminatedDoubleStranded

#
# other files
#
ref_genome := /mnt/solexa/Genomes/hg19_evan/whole_genome.fa

lippold_regions_bed := $(input_dir)/lippold_regions.bed
exome_regions_bed   := $(input_dir)/exome_regions.bed
lippold_sites_bed   := $(input_dir)/lippold_sites.bed
exome_sites_bed     := $(input_dir)/exome_sites.bed


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


init: $(data_dirs) $(bteam_info)

bams: $(data_dirs) $(lippold_bams) $(lippold_bais) $(exome_bams) $(exome_bais) $(denisova_bams)

genotypes: $(data_dirs) $(comb_merged_vcf) $(comb_merged_vcf).tbi $(lippold_merged_vcf) $(lippold_merged_vcf).tbi $(exome_merged_vcf) $(exome_merged_vcf).tbi $(dpfilt_comb_merged_vcf) $(dpfilt_comb_merged_vcf).tbi $(dpfilt_lippold_merged_vcf) $(dpfilt_lippold_merged_vcf).tbi $(dpfilt_exome_merged_vcf) $(dpfilt_exome_merged_vcf).tbi

alignments: $(data_dirs) $(lippold_fastas) $(exome_fastas) $(combined_fastas)

ancient_features: $(data_dirs)
	jupyter nbconvert $(nb_ancient_features) --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $(nb_ancient_features)

coverage_analysis: $(data_dirs)
	jupyter nbconvert $(nb_coverage_analysis) --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $(nb_coverage_analysis)

damage_patterns: $(bam_dir)/lippold_sidron.bam $(bam_dir)/exome_sidron.bam $(tmp_dir)/sidron_rmdup_len35mapq37_sorted.bam
	@cd $(figures_dir); \
	for bam in $^; do \
		python3 $(bam_plotdamage) --bam ../$$bam; \
	done


#
# BAM processing
#
$(bam_dir)/%_mez2.bam: $(input_dir)/%_regions.bed
	bedtools intersect -a /mnt/expressions/mateja/Late_Neandertals/Final_complete_dataset/Merged_per_individual_L35MQ0/Mezmaiskaya2_final.bam -b $< -sorted \
		> $@

$(bam_dir)/%_spy.bam: $(input_dir)/%_regions.bed
	bedtools intersect -a /mnt/expressions/mateja/Late_Neandertals/Final_complete_dataset/Merged_per_individual_L35MQ0/Spy_final.bam -b $< -sorted \
		> $@

$(bam_dir)/lippold_sidron.bam: $(lippold_regions_bed)
	jupyter nbconvert $(nb_sidron_processing) --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $(nb_sidron_processing); \
	mv $@ $@_tmp; \
	$(decrease_bquals) -n 5 $@_tmp $@; \
	rm $@_tmp

$(bam_dir)/exome_sidron.bam: $(exome_regions_bed) $(tmp_dir)/whole_exome.bam
	bedtools intersect -a $(tmp_dir)/whole_exome.bam -b $< -sorted \
		> $@

$(bam_dir)/%_a00_1.bam: $(input_dir)/%_regions.bed $(tmp_dir)/GRC13292545.chrY.bam
	bedtools intersect -a $(tmp_dir)/GRC13292545.chrY.bam -b $< \
		> $@; \
	samtools index $@

$(bam_dir)/%_a00_2.bam: $(input_dir)/%_regions.bed $(tmp_dir)/GRC13292546.chrY.bam
	bedtools intersect -a $(tmp_dir)/GRC13292546.chrY.bam -b $< \
		> $@; \
	samtools index $@

$(tmp_dir)/whole_exome.bam:
	curl http://cdna.eva.mpg.de/neandertal/exomes/BAM/Sidron_exome_hg19_1000g_LowQualDeamination.md.bam -o $@; \
	samtools index $@

$(tmp_dir)/GRC13292545.chrY.bam:
	cd $(tmp_dir); curl -O http://evolbio.ut.ee/chrY/GRC13292545.chrY.bam

$(tmp_dir)/GRC13292546.chrY.bam:
	cd $(tmp_dir); curl -O http://evolbio.ut.ee/chrY/GRC13292546.chrY.bam


$(denisova_bams): $(nb_denisova_processing) $(lippold_regions_bed)
	jupyter nbconvert $< --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $<; \touch $(denisova_bams)


#
# VCF processing
#
$(vcf_dir)/lippold_chimp.vcf.gz $(vcf_dir)/exome_chimp.vcf.gz: $(lippold_sites_bed) $(exome_sites_bed)
	jupyter nbconvert $(nb_chimpanzee_genotypes) --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $(nb_chimpanzee_genotypes); \
	touch $@

$(vcf_dir)/%_sidron.vcf.gz: $(input_dir)/%_regions.bed $(bam_dir)/%_sidron.bam
	samtools mpileup -l $(word 1,$^) -t DP -A -Q 20 -q 30 -u -f $(ref_genome) $(word 2,$^) \
		| bcftools call --ploidy 1 -m -V indels -Oz \
		| bcftools reheader -s <(echo -e "ElSidron"| cat) -o $@

$(vcf_dir)/%_sidron_cons.vcf.gz: $(bam_dir)/%_sidron.bam $(input_dir)/%_regions.bed
	python3 $(bam_sample) \
		--bam $(word 1,$^) --bed $(word 2,$^) \
		--ref $(ref_genome) --format VCF --sample-name ElSidronCons \
		--method consensus --minbq 20 --mincov 3 \
	| bgzip \
	> $@

$(vcf_dir)/%_qiaomei_sidron.vcf.gz: $(input_dir)/%_regions.bed
	 bcftools view -V indels /mnt/454/Carbon_beast_QM/Y_Sidron_TY/1_Extended_VCF/Sidron.hg19_evan.Y.mod.vcf.gz -R $< \
		| grep -v "LowQual" \
		| sed 's/0\/0/0/; s/1\/1/1/; s/\.\/\./\./' \
		| grep -v "0\/1" \
		| bcftools annotate -x INFO,FORMAT -Oz -o $@

$(vcf_dir)/%_mez2.vcf.gz: $(bam_dir)/%_mez2.bam $(input_dir)/%_regions.bed
	python3 $(bam_sample) \
		--bam $(word 1,$^) --bed $(word 2,$^) \
		--ref $(ref_genome) --format VCF --sample-name Mez2 \
		--strand-check non-USER_all --method majority \
	| bgzip \
	> $@

$(vcf_dir)/%_spy.vcf.gz: $(bam_dir)/%_spy.bam $(input_dir)/lippold_regions.bed
	python3 $(bam_sample) \
		--bam $(word 1,$^) --bed $(word 2,$^) \
		--ref $(ref_genome) --format VCF --sample-name Spy \
		--strand-check non-USER_all --method majority \
	| bgzip \
	> $@

$(vcf_dir)/%_a00_1.vcf.gz: $(input_dir)/%_regions.bed $(bam_dir)/%_a00_1.bam
	samtools mpileup -l $(word 1,$^) -t DP -A -Q 20 -q 30 -u -f $(ref_genome) $(word 2,$^) \
		| bcftools call --ploidy 1 -m -V indels -Oz \
		| bcftools reheader -s <(echo -e "A00_1"| cat) -o $@

$(vcf_dir)/%_a00_2.vcf.gz: $(input_dir)/%_regions.bed $(bam_dir)/%_a00_2.bam
	samtools mpileup -l $(word 1,$^) -t DP -A -Q 20 -q 30 -u -f $(ref_genome) $(word 2,$^) \
		| bcftools call --ploidy 1 -m -V indels -Oz \
		| bcftools reheader -s <(echo -e "A00_2"| cat) -o $@

$(vcf_dir)/lippold_bteam_%.vcf.gz: $(bteam_info) $(lippold_regions_bed)
	id=$(subst .vcf.gz,,$(subst lippold_bteam_,,$(notdir $@))); \
	$(process_bteam_vcf) $$id $^ $@

$(vcf_dir)/exome_bteam_%.vcf.gz: $(bteam_info) $(exome_regions_bed)
	id=$(subst .vcf.gz,,$(subst exome_bteam_,,$(notdir $@))); \
	$(process_bteam_vcf) $$id $^ $@


#
# merged VCFs
#
$(vcf_dir)/merged_lippold.vcf.gz: $(lippold_vcfs) $(lippold_tbis)
	bcftools merge -m all $(lippold_vcfs) \
		| bcftools view -M2 \
		| bcftools annotate -x INFO -Oz -o $@

$(vcf_dir)/merged_exome.vcf.gz: $(exome_vcfs) $(exome_tbis)
	bcftools merge -m all $(exome_vcfs) \
		| bcftools view -M2 \
		| bcftools annotate -x INFO -Oz -o $@

$(comb_merged_vcf): $(exome_merged_vcf) $(lippold_merged_vcf)
	bcftools concat $^ \
		| bcftools annotate -x INFO -Oz -o $@_unsorted; \
	zgrep '^#' $@_unsorted > $@_tmp; \
	zgrep -v '^#' $@_unsorted \
		| sort -k1,1n -k2,2n \
		>> $@_tmp; \
	bgzip -c $@_tmp > $@; \
	rm $@_unsorted $@_tmp

$(vcf_dir)/dpfilt_%.vcf.gz: $(vcf_dir)/%.vcf.gz
	vcftools --gzvcf $< --out $@_tmp --minDP 3 --recode; \
	bgzip $@_tmp.recode.vcf; \
	bcftools view -s ^Chimp,Mez2,Spy $@_tmp.recode.vcf.gz -Oz -o $@; \
	rm $@_tmp.recode.vcf.gz $@_tmp.log

#
# index files
#
$(bam_dir)/%.bai: $(bam_dir)/%
	samtools index $<

$(vcf_dir)/%.vcf.gz.tbi: $(vcf_dir)/%.vcf.gz
	tabix -f $<


#
# FASTA generation
#
sample_ids := ElSidron A00_1 A00_2 $(bteam_names)

$(fasta_dir)/%_chimp_nea_bteam.fa: $(vcf_dir)/merged_%.vcf.gz
	python $(src_dir)/vcf_to_fasta.py --vcf-file $< --fasta-file $@ --chrom Y --sample-names Chimp Mez2 Spy $(sample_ids)

$(fasta_dir)/%_sidron_bteam.fa: $(vcf_dir)/merged_%.vcf.gz
	python $(src_dir)/vcf_to_fasta.py --vcf-file $< --fasta-file $@ --chrom Y --sample-names $(sample_ids)

$(fasta_dir)/%_bteam.fa: $(vcf_dir)/merged_%.vcf.gz
	python $(src_dir)/vcf_to_fasta.py --vcf-file $< --fasta-file $@ --chrom Y --sample-names A00_1 A00_2 $(bteam_names)


#
# other things
#
$(lippold_regions_bed):
	bedtools intersect -a /mnt/454/Carbon_beast_QM/QF_chrY_region.bed \
		-b /mnt/454/HCNDCAM/Hengs_Alignability_Filter/hs37m_filt35_50.bed.gz \
		| sort -k1,1n -k2,2n \
		> $@

$(exome_regions_bed):
	curl http://www.cell.com/cms/attachment/2052899616/2060015784/mmc2.zip \
		| gunzip -c \
		> $@

$(input_dir)/%_sites.bed: $(input_dir)/%_regions.bed
	python $(src_dir)/sites_in_bed.py --bed-file $< --output-file $@ --format BED

$(bteam_info):
	grep '^SS.*M$$' /mnt/454/HighCovNeandertalGenome/1_Extended_VCF/Individuals.txt \
		| tr '\t' ' ' \
		| tr -s ' ' \
		| tr ' ' '\t' \
		| cut -f1,2,3 \
		> $@

$(input_dir)/coords_of_Y_exons.bed:
	curl http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz \
		| gunzip -c \
		| grep 'chrY' \
		| awk -vOFS='\t' '{ if ($$6 != $$7) { print $$2, $$6, $$7 }}' \
		| sed 's/chr//' \
		| sort -k1,1n -k2,2n \
		| bedtools merge -i stdin \
		> $@

$(data_dirs):
	mkdir -p $@


clean:
	rm -rf $(data_dirs)
