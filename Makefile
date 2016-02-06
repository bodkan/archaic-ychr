doc_dir := doc

bam_dir := bam
input_dir := input
figures_dir := figures
output_dir := output
tmp_dir := tmp
data_dirs := $(bam_dir) $(figures_dir) $(input_dir) $(output_dir) $(tmp_dir)

targets_bed := $(input_dir)/target_regions.bed

sidron_bam      := $(bam_dir)/sidron_ontarget.bam
den8_bam        := $(bam_dir)/den8_ontarget.bam
deam_den8_bam   := $(bam_dir)/deam_den8_ontarget.bam

nb_sidron_processing := $(doc_dir)/processing_of_El_Sidron_data.ipynb
nb_den8_processing := $(doc_dir)/processing_of_Denisova_8_data.ipynb
nb_adna_features := $(doc_dir)/analyze_aDNA_features.ipynb



process_bams: $(data_dirs) $(sidron_bam) $(den8_bam) $(deam_den8_bam)

analyze_bams:
	jupyter nbconvert $(sidron_processing_notebook) --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $(sidron_processing_notebook)



$(sidron_bam): $(targets_bed)
	jupyter nbconvert $(nb_sidron_processing) --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $(nb_sidron_processing)

$(den8_bam) $(den8_deam_bam): $(targets_bed)
	jupyter nbconvert $(nb_den8_processing) --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $(nb_den8_processing)



$(targets_bed):
	cp /mnt/454/Carbon_beast_QM/QF_chrY_region.bed $@

$(data_dirs):
	mkdir -p $@



clean:
	rm -rf $(figures_dir) $(input_dir) $(output_dir)

clean_all:
	rm -rf $(data_dirs)
