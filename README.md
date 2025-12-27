### Activate conda env for qiime


	conda activate qiime2-amplicon-xxxx.x



### For qiime2 demultiplexed paired ends reads, binned quality scores

Create manifest.txt file (tab delimited), add some columns for metadata on sample if possible, e.g. genotype etc**

**1st column - sample-id**


	ls *1.fq.gz | cut -c x-x > name.txt


**2nd column - forward-absolute-filepath**


	realpath *1.fq.gz > forward.txt


**3rd column - reverse-absolute-filepath**


	realpath *2.fq.gz > reverse.txt


**For other columns, add any other relevant metadata such as batch/cage etc.**

**Info on formating metadata here: https://docs.qiime2.org/2024.10/tutorials/metadata/**


### Import our reads

	mkdir output

	qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifest.txt --output-path output/demux --input-format PairedEndFastqManifestPhred33V2


### Visualize sequence quality with fastQC/multiQC. qiime2 doesnt work for binned QC scores

	mkdir fastqc

	fastqc -o ./fastqc *.fq.gz

	multiqc ./fastqc

 
### Trim sequences off the primer sequence (below is based on v3v4 - Azenta primer seq)

	cd output

	qiime cutadapt trim-paired --i-demultiplexed-sequences demux.qza --p-front-f ACTCCTACGGGAGGCAGCAG --p-front-r GGACTACHVGGGTWTCTAAT --o-trimmed-sequences trimmed-demux.qza --p-match-read-wildcards --p-discard-untrimmed --verbose

 
### Visualize the trimmed sequence quality with fastQC/multiQC
   
	qiime tools extract --input-path trimmed-demux.qza --output-path trimmed_qc/
 	qiime demux summarize --i-data trimmed-demux.qza --o-visualization trimmed-demux.qzv

	repeat 3. but change directory


### DADA2 (trunc len depends on the QC done in step 5, use --p-trunc-len-f/r to trim forward/reverse read based on quality scores if needed)
### *add --p-n-threads x if PC can support the ram/cpu usage.

	qiime dada2 denoise-paired --i-demultiplexed-seqs trimmed-demux.qza --p-trunc-len-f 0 --p-trunc-len-r 0 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza --o-base-transition-stats base-transition-stats.qza --verbose --p-n-threads 6 --p-no-hashed-feature-ids


### View summary of statistics

	qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file <path>/manifest.txt

	qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv

	qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv


### Generate a tree for phylogenetic diversity analyses

	qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.qza --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza


### Construct reference database based on silva 16S rRNA gene database, that is amplicon/primer specific
### Refer to (https://forum.qiime2.org/t/processing-filtering-and-evaluating-the-silva-database-and-other-reference-sequence-data-with-rescript/15494)
### For silva based taxonomy

	mkdir silva
 
	cd silva
 
	wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/taxonomy/tax_slv_ssu_138.2.txt.gz
 
	wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.2.txt.gz
 
	wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/taxonomy/tax_slv_ssu_138.2.tre.gz
 
	wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz
 
	qiime tools import --type 'FeatureData[SILVATaxonomy]' --input-path tax_slv_ssu_138.2.txt --output-path taxranks-silva-138.2-ssu-nr99.qza
 
	qiime tools import --type 'FeatureData[SILVATaxidMap]' --input-path taxmap_slv_ssu_ref_nr_138.2.txt --output-path taxmap-silva-138.2-ssu-nr99.qza
 
	qiime tools import --type 'Phylogeny[Rooted]' --input-path tax_slv_ssu_138.2.tre --output-path taxtree-silva-138.2-nr99.qza
 
	qiime tools import --type 'FeatureData[RNASequence]' --input-path SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta --output-path silva-138.2-ssu-nr99-rna-seqs.qza
 
	qiime rescript reverse-transcribe --i-rna-sequences silva-138.2-ssu-nr99-rna-seqs.qza --o-dna-sequences silva-138.2-ssu-nr99-seqs.qza
 
	qiime rescript parse-silva-taxonomy --i-taxonomy-tree taxtree-silva-138.2-nr99.qza --i-taxonomy-map taxmap-silva-138.2-ssu-nr99.qza --i-taxonomy-ranks taxranks-silva-138.2-ssu-nr99.qza --o-taxonomy silva-138.2-ssu-nr99-tax.qza --p-rank-propagation
 
	qiime rescript cull-seqs --i-sequences silva-138.2-ssu-nr99-seqs.qza --o-clean-sequences silva-138.2-ssu-nr99-seqs-cleaned.qza --p-n-jobs 72
 
	qiime rescript filter-seqs-length-by-taxon --i-sequences silva-138.2-ssu-nr99-seqs-cleaned.qza --i-taxonomy silva-138.2-ssu-nr99-tax.qza --p-labels Archaea Bacteria Eukaryota --p-min-lens 900 1200 1400 --o-filtered-seqs silva-138.2-ssu-nr99-seqs-filt.qza --o-discarded-seqs silva-138.2-ssu-nr99-seqs-discard.qza
 
	qiime rescript dereplicate --i-sequences silva-138.2-ssu-nr99-seqs-filt.qza --i-taxa silva-138.2-ssu-nr99-tax.qza --p-mode 'uniq' --o-dereplicated-sequences silva-138.2-ssu-nr99-seqs-derep-uniq.qza --o-dereplicated-taxa silva-138.2-ssu-nr99-tax-derep-uniq.qza
 
	qiime feature-classifier extract-reads --i-sequences silva-138.2-ssu-nr99-seqs-derep-uniq.qza --p-f-primer ACTCCTACGGGAGGCAGCAG --p-r-primer GGACTACHVGGGTWTCTAAT --p-read-orientation 'forward' --o-reads silva-138.2-ssu-nr99-seqs-338f-806r.qza --p-n-jobs 72
 
	qiime rescript dereplicate --i-sequences silva-138.2-ssu-nr99-seqs-338f-806r.qza --i-taxa silva-138.2-ssu-nr99-tax-derep-uniq.qza --p-mode 'uniq' --o-dereplicated-sequences silva-138.2-ssu-nr99-seqs-338f-806r-uniq.qza --o-dereplicated-taxa  silva-138.2-ssu-nr99-tax-338f-806r-derep-uniq.qza
 
	qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads silva-138.2-ssu-nr99-seqs-338f-806r-uniq.qza --i-reference-taxonomy silva-138.2-ssu-nr99-tax-338f-806r-derep-uniq.qza --o-classifier silva-138.2-ssu-nr99-338f-806r-classifier.qza
 
	qiime feature-classifier classify-sklearn --i-classifier silva-138.2-ssu-nr99-338f-806r-classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza


 ### For greengenes2 based taxonomy
 ### Note: A new phylo tree has to be made based on gg2 ref-seq table
 
	mkdir greengenes2
 
	cd greengenes2
 
	wget https://ftp.microbio.me/greengenes_release/2024.09/2024.09.backbone.full-length.fna.qza
 
	wget https://ftp.microbio.me/greengenes_release/current/2024.09.backbone.tax.qza
 
	wget https://ftp.microbio.me/greengenes_release/current/2024.09.taxonomy.asv.nwk.qza
 
	qiime greengenes2 non-v4-16s --i-table <path>/table.qza --i-sequences <path>/rep-seqs.qza --i-backbone 2024.09.backbone.full-length.fna.qza --o-mapped-table gg2-filtered-table.qza --o-representatives gg2-filtered-rep-seqs.qza
 
	qiime greengenes2 taxonomy-from-table --i-reference-taxonomy 2024.09.taxonomy.asv.nwk.qza --i-table gg2-filtered-table.qza --o-classification taxonomy.qza

 	qiime phylogeny align-to-tree-mafft-fasttree --i-sequences gg2-filtered-rep-seqs.qza --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

### Visualize taxonomy
    
	qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv

    
### Proceed to R

### For picrust2

	qiime tools export --input-path rep-seqs.qza --output-path rep-seqs.fasta

	qiime tools export --input-path table.qza --output-path table.biom
 
	conda activate picrust2

 	picrust2_pipeline.py -s dna-sequences.fasta -i feature-table.biom -o picrust2_out_pipeline -p 1

 	qiime picrust2 full-pipeline --i-table table.qza --i-seq rep-seqs.qza --output-dir q2-picrust2_output --p-placement-tool epa-n --p-threads 1 --p-hsp-method m --p-max-nsti 2 --verbose

### Files required are: table.qza (or gg2 specific table), rooted-tree.qza, taxonomy.qza, manifest.txt
