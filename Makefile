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
all_bams := mez2.bam spy.bam elsidron.bam ustishim.bam a00.bam
exome_bams := $(addprefix $(bam_dir)/, $(addprefix exome_,$(all_bams)))

# VCF files
all_vcfs := mez2.vcf.gz spy.vcf.gz elsidron.vcf.gz ustishim.vcf.gz a00.vcf.gz
exome_vcfs := $(addprefix $(vcf_dir)/, $(addprefix exome_,$(all_vcfs)))

# FASTA files
fastas := chimp_nea_bteam.fa sidron_bteam.fa bteam.fa
exome_fastas := $(addprefix $(fasta_dir)/exome_,$(fastas))

# scripts
bam_sample := $(dep_dir)/bam-sample/bam-sample
run_nb := $(src_dir)/run_nb.sh

# coordinates
full_coord := $(coord_dir)/capture_full.bed
lippold_coord := $(coord_dir)/capture_lippold.bed
exome_coord := $(coord_dir)/capture_exome.bed
annot_coord := $(coord_dir)/cds.bed $(coord_dir)/phastcons.bed $(coord_dir)/genes.bed $(coord_dir)/pseudogenes.bed

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

init: $(dirs) $(full_coord) $(lippold_coord) $(exome_coord) $(annot_coord)

bam: $(dirs) $(full_bams) $(lippold_bams) $(exome_bams)

vcf: $(dirs) $(full_vcfs) $(lippold_vcfs) $(exome_vcfs)

fasta: $(dirs) $(full_fastas) $(lippold_fastas) $(exome_fastas)

damage_patterns:



#
# BAM processing
#
$(bam_dir)/%_mez2.bam: $(input_dir)/%_capture.bed
	bedtools intersect -a /mnt/expressions/mateja/Late_Neandertals/Final_complete_dataset/Merged_per_individual_L35MQ0/Mezmaiskaya2_final.bam -b $< -sorted \
		> $@

$(bam_dir)/%_spy.bam: $(input_dir)/%_capture.bed
	bedtools intersect -a /mnt/expressions/mateja/Late_Neandertals/Final_complete_dataset/Merged_per_individual_L35MQ0/Spy_final.bam -b $< -sorted \
		> $@

$(bam_dir)/lippold_sidron.bam: $(lippold_coord)
	jupyter nbconvert $(nb_sidron_processing) --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $(nb_sidron_processing); \
	mv $@ $@_tmp; \
	$(decrease_bquals) -n 5 $@_tmp $@; \
	rm $@_tmp

$(bam_dir)/exome_sidron.bam: $(exome_coord) $(tmp_dir)/whole_exome.bam
	bedtools intersect -a $(tmp_dir)/whole_exome.bam -b $< -sorted \
		> $@_unfilt; \
	bam-rmdup -l 35 -q 37 -r -o $@ $@_unfilt

$(bam_dir)/%_ust_ishim.bam: $(input_dir)/%_capture.bed
	bedtools intersect -a /mnt/expressions/mp/y-selection/bam/y_ustishim.bam -b $< \
		> $@

$(bam_dir)/%_a00.bam: $(bam_dir)/%_a00_1.bam $(bam_dir)/%_a00_2.bam
	samtools merge $@ $^

$(bam_dir)/%_a00_1.bam: $(input_dir)/%_capture.bed $(tmp_dir)/GRC13292545.chrY.bam
	bedtools intersect -a $(tmp_dir)/GRC13292545.chrY.bam -b $< \
		> $@; \
	samtools index $@

$(bam_dir)/%_a00_2.bam: $(input_dir)/%_capture.bed $(tmp_dir)/GRC13292546.chrY.bam
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



#
# VCF processing
#
$(vcf_dir)/lippold_chimp.vcf.gz $(vcf_dir)/exome_chimp.vcf.gz: $(lippold_sites_bed) $(exome_sites_bed)
	jupyter nbconvert $(nb_chimpanzee_genotypes) --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $(nb_chimpanzee_genotypes); \
	touch $@

$(vcf_dir)/%_sidron.vcf.gz: $(input_dir)/%_capture.bed $(bam_dir)/%_sidron.bam
	samtools mpileup -l $(word 1,$^) -t DP -A -Q 20 -q 30 -u -f $(ref_genome) $(word 2,$^) \
		| bcftools call --ploidy 1 -m -V indels -Oz \
		| bcftools reheader -s <(echo -e "ElSidron"| cat) -o $@

$(vcf_dir)/%_sidron_cons.vcf.gz: $(bam_dir)/%_sidron.bam $(input_dir)/%_capture.bed
	python3 $(bam_sample) \
		--bam $(word 1,$^) --bed $(word 2,$^) \
		--ref $(ref_genome) --format VCF --sample-name ElSidronCons \
		--method consensus --minbq 20 --mincov 3 \
	| bgzip \
	> $@

$(vcf_dir)/%_qiaomei_sidron.vcf.gz: $(input_dir)/%_capture.bed
	 bcftools view -V indels /mnt/454/Carbon_beast_QM/Y_Sidron_TY/1_Extended_VCF/Sidron.hg19_evan.Y.mod.vcf.gz -R $< \
		| grep -v "LowQual" \
		| sed 's/0\/0/0/; s/1\/1/1/; s/\.\/\./\./' \
		| grep -v "0\/1" \
		| bcftools annotate -x INFO,FORMAT -Oz -o $@

$(vcf_dir)/%_mez2.vcf.gz: $(bam_dir)/%_mez2.bam $(input_dir)/%_capture.bed
	python3 $(bam_sample) \
		--bam $(word 1,$^) --bed $(word 2,$^) \
		--ref $(ref_genome) --format VCF --sample-name Mez2 \
		--strand-check non-USER_all --method majority \
	| bgzip \
	> $@

$(vcf_dir)/%_spy.vcf.gz: $(bam_dir)/%_spy.bam $(input_dir)/lippold_capture.bed
	python3 $(bam_sample) \
		--bam $(word 1,$^) --bed $(word 2,$^) \
		--ref $(ref_genome) --format VCF --sample-name Spy \
		--strand-check non-USER_all --method majority \
	| bgzip \
	> $@

$(vcf_dir)/%_ust_ishim.vcf.gz: $(input_dir)/%_capture.bed $(bam_dir)/%_ust_ishim.bam
	samtools mpileup -l $(word 1,$^) -t DP -A -Q 20 -q 30 -u -f $(ref_genome) $(word 2,$^) \
		| bcftools call --ploidy 1 -m -V indels -Oz -o $@

$(vcf_dir)/%_a00.vcf.gz: $(input_dir)/%_capture.bed $(bam_dir)/%_a00.bam
	samtools mpileup -l $(word 1,$^) -t DP -A -Q 20 -q 30 -u -f $(ref_genome) $(word 2,$^) \
		| bcftools call --ploidy 1 -m -V indels -Oz \
		| bcftools reheader -s <(echo -e "A00"| cat) -o $@

$(vcf_dir)/lippold_bteam_%.vcf.gz: $(bteam_info) $(lippold_coord)
	id=$(subst .vcf.gz,,$(subst lippold_bteam_,,$(notdir $@))); \
	$(process_bteam_vcf) $$id $^ $@

$(vcf_dir)/exome_bteam_%.vcf.gz: $(bteam_info) $(exome_coord)
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
# FASTA alignments for three different capture sets
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
$(lippold_coord):
	scp bionc11.eva.mpg.de:/mnt/genotyping/sendru/basti_design.bed $@

# Y chromosome capture regions designed by Qiaomei (~6Mb)
$(full_coord):
	scp bionc11.eva.mpg.de:/mnt/454/Carbon_beast_QM/array_2015_0729/array_order/Y.filt35_50_SRepeat_100.bed $(tmp_dir)/; \
	perl -lane 'print $$F[0] . "\t" . $$F[1] . "\t" . $$F[2]' $(tmp_dir)/Y.filt35_50_SRepeat_100.bed \
	    > $@

# Y chromosome exome capture regions
$(exome_coord):
	wget http://www.cell.com/cms/attachment/2052899616/2060015784/mmc2.zip
	unzip mmc2.zip; rm mmc2.zip
	mv ajhg2064mmc2_V1.txt $@

# functional annotation coordinates
$(coord_dir)/cds.bed $(coord_dir)/phastcons.bed $(coord_dir)/genes.bed $(coord_dir)/pseudogenes.bed:
	$(run_nb) notebooks/annotated_regions.ipynb



#
# other dependencies
#
$(bam_sample):
	git clone https://github.com/bodkan/bam-sample $(dep_dir)/bam-sample

$(dirs):
	mkdir -p $@

clean:
	rm -rf $(dirs)
