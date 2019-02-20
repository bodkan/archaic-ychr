SHELL := /bin/bash

# directories
data_dir := data
bam_dir := $(data_dir)/bam
fasta_dir := $(data_dir)/fasta
vcf_dir := $(data_dir)/vcf
coord_dir := $(data_dir)/coord
tmp_dir := tmp
fig_dir := fig
src_dir := src
dirs := $(data_dir) $(bam_dir) $(vcf_dir) $(fasta_dir) $(coord_dir) $(fig_dir) $(tmp_dir)

# BAM files
sgdp_bams :=  S_French-1.bam S_Sardinian-1.bam S_Han-2.bam S_Dai-2.bam S_Papuan-2.bam S_Karitiana-1.bam S_Dinka-1.bam S_Mbuti-1.bam S_Yoruba-2.bam S_Mandenka-1.bam
published_bams := $(sgdp_bams) ustishim.bam a00_1.bam a00_2.bam kk1.bam mota.bam bichon.bam loschbour.bam
full_bams := $(addprefix $(bam_dir)/, $(addprefix full_, spy1.bam mez2.bam denisova8.bam $(published_bams)))
lippold_bams := $(addprefix $(bam_dir)/, $(addprefix lippold_, elsidron2.bam $(published_bams)))
exome_bams := $(addprefix $(bam_dir)/, $(addprefix exome_, elsidron1.bam $(published_bams)))

# VCF files
published_vcfs := $(subst .bam,.vcf.gz, $(published_bams))
full_vcfs := $(addprefix $(vcf_dir)/, $(addprefix full_, spy1.vcf.gz mez2.vcf.gz denisova8.vcf.gz $(published_vcfs)))
lippold_vcfs := $(addprefix $(vcf_dir)/, $(addprefix lippold_, elsidron2.vcf.gz $(published_vcfs)))
exome_vcfs := $(addprefix $(vcf_dir)/, $(addprefix exome_, elsidron1.vcf.gz $(published_vcfs)))

full_vcf := $(vcf_dir)/merged_full.vcf.gz
lippold_vcf := $(vcf_dir)/merged_lippold.vcf.gz
exome_vcf := $(vcf_dir)/merged_exome.vcf.gz

# FASTA files
fastas := chimp_nea_bteam.fa sidron_bteam.fa bteam.fa
exome_fastas := $(addprefix $(fasta_dir)/exome_,$(fastas))

# scripts
bam_sample := /mnt/expressions/mp/bam-sample/bam-sample
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

init: $(dirs) $(full_bed) $(lippold_bed) $(exome_bed) $(annot_bed)

bam: $(dirs) $(full_bams) $(lippold_bams) $(exome_bams)

vcf: $(dirs) $(full_vcf) $(lippold_vcf) $(exome_vcf)

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
	cp /mnt/genotyping/sk_pipelines/datasets/Mallick2016_SGDP_Ychromosome/$(basename $(notdir $@)).Y.bam $@.all
	bedtools intersect -a $@.all -b $(map_filter) > $@; rm $@.all

# A00 Y
$(tmp_dir)/a00_1.bam:
	cd $(tmp_dir); curl -O http://evolbio.ut.ee/chrY/GRC13292545.chrY.bam
	bedtools intersect -a $(tmp_dir)/GRC13292545.chrY.bam -b $(map_filter) > $@

$(tmp_dir)/a00_2.bam:
	cd $(tmp_dir); curl -O http://evolbio.ut.ee/chrY/GRC13292546.chrY.bam
	bedtools intersect -a $(tmp_dir)/GRC13292546.chrY.bam -b $(map_filter) > $@

$(tmp_dir)/elsidron1.bam:
	cd $(tmp_dir); curl -O http://cdna.eva.mpg.de/neandertal/exomes/BAM/Sidron_exome_hg19_1000g_LowQualDeamination.md.bam; \
	    $(analyze_bam) -qual 25 -minlength 35 Sidron_exome_hg19_1000g_LowQualDeamination.md.bam
	bedtools intersect -a $(tmp_dir)/Sidron_exome_hg19_1000g_LowQualDeamination.md.uniq.L35MQ25.bam -b $(map_filter) > $@
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
	samtools view -h -b /mnt/expressions/mp/Archive/y-selection/tmp/KK1_sort_rmdup_merge_IR_q30_mapDamage.bam Y -o $(tmp_dir)/KK1.Y.bam
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 KK1.Y.bam
	bedtools intersect -a $(tmp_dir)/KK1.Y.uniq.L35MQ25.bam -b $(map_filter) > $@
	samtools index $@

$(tmp_dir)/bichon.bam:
	samtools view -h -b /mnt/expressions/mp/Archive/y-selection/tmp/Bichon.sort.rmdup.IR.q30.mapDamage.bam Y -o $(tmp_dir)/Bichon.Y.bam
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 Bichon.Y.bam
	bedtools intersect -a $(tmp_dir)/Bichon.Y.uniq.L35MQ25.bam -b $(map_filter) > $@
	samtools index $@

$(tmp_dir)/mota.bam:
	samtools view -h -b /mnt/expressions/mp/Archive/y-selection/tmp/GB20_sort_merge_dedup_l30_IR_q30_mapDamage.bam Y -o $(tmp_dir)/Mota.Y.bam
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 Mota.Y.bam
	bedtools intersect -a $(tmp_dir)/Mota.Y.uniq.L35MQ25.bam -b $(map_filter) > $@
	samtools index $@

$(tmp_dir)/loschbour.bam:
	samtools view -h -b /mnt/expressions/mp/Archive/y-selection/tmp/Loschbour.hg19_1000g.bam Y -o $(tmp_dir)/Loschbour.Y.bam
	cd $(tmp_dir); $(analyze_bam) -qual 25 -minlength 35 Loschbour.Y.bam
	bedtools intersect -a $(tmp_dir)/Loschbour.Y.uniq.L35MQ25.bam -b $(map_filter) > $@
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

$(vcf_dir)/merged_full.vcf.gz: $(vcf_dir)/full_chimp.vcf.gz $(full_vcfs)
	bcftools merge $^ | bcftools annotate -x INFO -Oz -o $@

$(vcf_dir)/merged_lippold.vcf.gz: $(vcf_dir)/lippold_chimp.vcf.gz $(lippold_vcfs)
	bcftools merge $^ | bcftools annotate -x INFO -Oz -o $@

$(vcf_dir)/merged_exome.vcf.gz: $(vcf_dir)/exome_chimp.vcf.gz $(exome_vcfs)
	bcftools merge $^ | bcftools annotate -x INFO -Oz -o $@

$(vcf_dir)/%_chimp.vcf.gz: $(coord_dir)/capture_%.pos
	$(src_dir)/chimp_vcf.sh $< $(basename $@)
	bgzip $(basename $@)
	tabix $@

$(vcf_dir)/%.vcf.gz: $(bam_dir)/%.bam
	$(bam_sample) --bam $< --ref $(ref_genome) --strategy consensus --format vcf \
	    --sample-name $(shell echo $(basename $(notdir $<)) | sed 's/^[a-z]*_//') --output $(basename $(basename $@))
	bgzip $(basename $@)
	tabix $@



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
$(coord_dir)/capture_%.pos: $(coord_dir)/capture_%.bed
	$(src_dir)/sites_in_bed.py --bed $< --format pos --output $@



$(dirs):
	mkdir -p $@

clean:
	rm -rf $(dirs)
