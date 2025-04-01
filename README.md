### Activate conda env for qiime

``
	conda activate qiime2-amplicon-xxxx.x
``


### For qiime2 demultiplexed paired ends reads, binned quality scores

Create manifest.txt file (tab delimited), add some columns for metadata on sample if possible, e.g. genotype etc**

**1st column - sample-id**

``
	ls *1.fq.gz | cut -c x-x > name.txt
``

**2nd column - forward-absolute-filepath**

``
	realpath *1.fq.gz > forward.txt
``

**3rd column - reverse-absolute-filepath**

``
	realpath *2.fq.gz > reverse.txt
``

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

	qiime cutadapt trim-paired --i-demultiplexed-sequences demux.qza --p-front-f ACTCCTACGGGAGGCAGCAG --p-front-r GGACTACHVGGGTWTCTAAT --o-trimmed-sequences trimmed-demux.qza --p-match-read-wildcards --verbose

 
### Visualize the trimmed sequence quality with fastQC/multiQC
   
	qiime tools extract --input-path trimmed-demux.qza --output-path trimmed_qc/

	repeat 3. but change directory


### dada2 (trunc len depends on the QC done in step 5, use --p-trunc-len-f/r for forward/reverse read) *add --p-n-threads 0 if PC can support the ram usage. Change the number accordingly to the sequence length.

	qiime dada2 denoise-paired --i-demultiplexed-seqs trimmed-demux.qza --p-trunc-len-f 0 --p-trunc-len-r 0 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza --verbose --p-n-threads 72 --p-no-hashed-feature-ids


### View summary of statistics

	qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file <path>/manifest.txt

	qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv

	qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv


### Generate a tree for phylogenetic diversity analyses

	qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.qza --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza


### Construct reference database based on silva 16S rRNA gene database, that is amplicon/primer specific - only need to be done if using different primer and a new database (https://forum.qiime2.org/t/processing-filtering-and-evaluating-the-silva-database-and-other-reference-sequence-data-with-rescript/15494) -> for silva based taxonomy

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
 
	qiime feature-classifier classify-sklearn --i-classifier silva-138.2-ssu-nr99-338f-806r-classifier.qza --i-reads <path>/rep-seqs.qza --o-classification taxonomy.qza


 ### For greengenes2 based taxonomy
 
	mkdir greengenes2
 
	cd greengenes2
 
	wget https://ftp.microbio.me/greengenes_release/2024.09/2024.09.backbone.full-length.fna.qza
 
	wget https://ftp.microbio.me/greengenes_release/current/2024.09.backbone.tax.qza
 
	wget https://ftp.microbio.me/greengenes_release/current/2024.09.taxonomy.asv.nwk.qza
 
	qiime greengenes2 non-v4-16s --i-table <path>/table.qza --i-sequences <path>/rep-seqs.qza --i-backbone 2024.09.backbone.full-length.fna.qza --o-mapped-table gg2-filtered-table.qza --o-representatives gg2-filtered-rep-seqs.qza
 
	qiime greengenes2 taxonomy-from-table --i-reference-taxonomy 2024.09.taxonomy.asv.nwk.qza --i-table gg2-filtered-table.qza --o-classification taxonomy.qza

### Visualize taxonomy
    
	qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv

    
### Proceed to R

Files required are: table.qza, rooted-tree.qza, taxonomy.qza, manifest.txt

	library(speedyseq)
	library(qiime2R)
	library(microViz)
	library(data.table)
	library(tidyverse)
	library(patchwork)
	library(ggrepel)
	library(ComplexHeatmap)
	library(phyloseq)
	library(ggnewscale)
	library(lefser)
	library(mia)
	library(ggforce)
	library(ape)
	library(doParallel)
	library(ConQuR)
	library(mirlyn)
 	
  	# Create dir to hold plots
	dir.create("./visualization")
	plots.dir <- "./visualization"


	# Store DNA Sequences
	dna <- Biostrings::DNAStringSet(taxa_names(ps))
	names(dna) <- taxa_names(ps)
	ps <- merge_phyloseq(ps, dna) 
	taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

	# Check for presence of batch effects
	batch.eff_wk0 <- ps %>%
  	tax_fix() %>%
  	ord_calc(method = "PCA") %>%
  	ord_plot(alpha = 0.6, size = 2, color = "Batch") +
  	theme_classic(12) +
  	geom_mark_ellipse(aes(color = Batch)) +
  	scale_y_reverse()
   
	batch.eff_wk0
 
	# Save plot
	ggsave(plot = batch.eff_wk0, path = plots.dir, filename = ("batch_eff.pdf"), 
       width = 14, height = 10, dpi = "retina")

	# Attempt to remove btch effects
	B <- as.data.frame(ps@otu_table)
	B <- t(B)
	B <- as.data.frame(B)

	# Extract batch ID from sample data
	batchid <- ps@sam_data$Batch
	summary(batchid)

	# Extract covariates
	D = ps@sam_data[, c('Genotype')] #covar
	summary(D)

	# Correct for batch effects using ConQuR package
	options(warn=-1) # required to call
	set.seed(2702)
	tuned_taxa <- Tune_ConQuR(tax_tab=B, batchid=batchid, covariates=D,
                            batch_ref_pool=c("1", "2"),
                            logistic_lasso_pool=F, 
                            quantile_type_pool=c("standard", "lasso"),
                            simple_match_pool=F,
                            lambda_quantile_pool=c(NA, "2p/n"),
                            interplt_pool=F,
                            frequencyL=0,
                            frequencyU=1,
                            num_core = 2)

       # Transpose the corrected matrix and convert it to a data frame
	taxa_correct2 <- t(tuned_taxa$tax_final)
	taxa_correct2 <- as.data.frame(taxa_correct2)

 	# Create new ASV table, taxonomy table, sample data, phy tree and refseq
	ASV <- otu_table(taxa_correct2, taxa_are_rows = TRUE)
	TAXA <- tax_table(ps)
	sampledata <- sample_data(ps)
	phytree <- phy_tree(ps)
	refseq <- refseq(ps)

	# repack the objects into a level 4 phyloseq structural data
	ps.bat <- phyloseq(ASV, TAXA, sampledata, phytree, refseq)

	 # Check post batch effect correction
	batch.eff.corr <- ps.bat %>%
  	tax_fix() %>%
  	ord_calc(method = "PCA") %>%
  	ord_plot(alpha = 0.6, size = 2, color = "Batch") +
	theme_classic(12) +
  	geom_mark_ellipse(aes(color = Batch)) +
  	scale_y_reverse()
   
	batch.eff.corr

	batch.compare <- batch.eff + batch.eff.corr

 	# Save plot
	ggsave(plot = batch.compare, path = plots.dir, filename = ("batch_compare.pdf"), 
       	width = 14, height = 5, dpi = "retina")

 	# Remove Unassigned in kingdom of tax table
	ps <- subset_taxa(ps, Kingdom != "Unassigned")
	ps.bat <- subset_taxa(ps.bat, Kingdom != "Unassigned")

 	# Remove Eukaroyta
	ps <- subset_taxa(ps, Kingdom != "d__Eukaryota")
	ps.bat <- subset_taxa(ps.bat, Kingdom != "d__Eukaryota")

 	# Remove NA in Phylum of tax table
	ps <- subset_taxa(ps, Phylum != "NA")
	ps.bat <- subset_taxa(ps.bat, Phylum != "NA")

 	# Fix NA values in the tax_table of other ranks by propagation
	ps <- tax_fix(ps)
	ps.bat <- tax_fix(ps.bat)

 	# Get read counts
	readcount_ps <- data.table(as(sample_data(ps.bat), "data.frame"),
                             TotalReads = sample_sums(ps.bat), 
                             keep.rownames = TRUE)
	write.csv(readcount_ps, "readcount_ps.csv")
	head(readcount_ps[order(readcount_ps$TotalReads), c("rn", "TotalReads")])

 	# Generate rarefaction curve to estimate the depth to rarefy
	ps <- orient_taxa(ps, as = "rows")
	rarefy_ps <- rarefy_whole_rep(ps, rep = 100, set.seed = 2702)
	rare_curve <- rarecurve(rarefy_ps, sample = "Sample")
	rare_curve
	ggsave(plot = rare_curve, path = plots.dir, filename = ("rarecurve.pdf"), 
       	width = 14, height = 10, dpi = "retina")

 	# Repeated rarefy
	mirl <- mirl((ps), rep = 1000, set.seed = 2702)
 	mirl_otu <- vector("list", length(mirl))
	for (i in 1:length(mirl)){
  	colnames(mirl[[i]]@otu_table) <- paste0(colnames(mirl[[i]]@otu_table), "_",i)
  	(mirl_otu[[i]] <- mirl[[i]]@otu_table)
	}
 
	mirl_rep_df <- do.call(cbind, mirl_otu)
	ps_id <- ps@sam_data@row.names
	average_counts <- vector("list", length(ps_id))
	for (i in 1:length(ps_id)){
  	sample_df <- mirl_rep_df[,grepl(ps_id[[i]], colnames(mirl_rep_df))]
  	sample_average <- data.frame(rowMeans(sample_df))
  	colnames(sample_average) <- ps_id[[i]]
 	average_counts[[i]] <- sample_average
	}
 
	average_count_df <- do.call(cbind, average_counts)
	ps.mirl <- ps
	average_count_dm <- data.matrix(average_count_df)
	ps.mirl@otu_table@.Data <- average_count_dm
