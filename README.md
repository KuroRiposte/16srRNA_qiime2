0. Activate conda env for qiime
	conda activate qiime2-amplicon-xxxx.x

For qiime2 demultiplexed paired ends reads, binned quality scores

1. Create manifest.txt file (tab delimited), add some columns for metadata on sample if possible, e.g. genotype etc
	1st column - sample-id
	ls *1.fq.gz | cut -c x-x > name.txt
	2nd column - forward-absolute-filepath
	realpath *1.fq.gz > forward.txt
	3rd column - reverse-absolute-filepath
	realpath *2.fq.gz > reverse.txt
	for other columns, add any other relevant metadata such as batch/cage etc. 
	info on formating metadata here: https://docs.qiime2.org/2024.10/tutorials/metadata/

2. Import our reads
	mkdir output
	qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifest.txt --output-path output/demux --input-format PairedEndFastqManifestPhred33V2

3. Visualize sequence quality with fastQC/multiQC. qiime2 doesnt work for binned QC scores
	mkdir fastqc
	fastqc -o ./fastqc *.fq.gz
	multiqc ./fastqc
	
4. Trim sequences off the primer sequence (below is based on v3v4 - Azenta primer seq)
	cd output
	qiime cutadapt trim-paired --i-demultiplexed-sequences demux.qza --p-front-f ACTCCTACGGGAGGCAGCAG --p-front-r GGACTACHVGGGTWTCTAAT --o-trimmed-sequences trimmed-demux.qza --p-match-read-wildcards --verbose
	
5. Visualize the trimmed sequence quality with fastQC/multiQC
	qiime tools extract --input-path trimmed-demux.qza --output-path trimmed_qc/
	blah blah repeat 3. but change directory

6. dada2 (trunc len depends on the QC done in step 5, use --p-trunc-len-f/r for forward/reverse read) *add --p-n-threads 0 if PC can support the ram usage. Change the number accordingly to the sequence length.
	qiime dada2 denoise-paired --i-demultiplexed-seqs trimmed-demux.qza --p-trunc-len-f 0 --p-trunc-len-r 225 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza --verbose --p-n-threads 6 --p-no-hashed-feature-ids

7. View summary of statistics
	qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file path to /manifest.txt
	qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv
	qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv
