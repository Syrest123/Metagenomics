# Qiime2 analysis of 16S data
![Qiime2 pipeline](https://murraycadzow.github.io/2021-obss-day4/fig/qiime2pipeline.png)
## Importing data into qiime
Importing data into qiime2 using a manifest file (sample-id\tforward-absolute-filepath\treverse-absolute-filepath)
```
#!/bin/env/bash
# Timing the entire pipeline
START="$( date +%s )"
```
```
# Importing data into qiime using the manifest format
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.tsv \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

```
Generates the summary of the distribution of the sequence quality at each position and also how many sequences were obtained per sample
```
qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization demux.qzv
```
## Quality control
```
qiime dada2 denoise-single \
  --i-demultiplexed-seqs paired-end-demux.qza \
  --p-trim-left 0 \
  --p-trunc-len 300 \
  --o-representative-sequences rep-seqs-dada2.qza \
  --o-table table-dada2.qza \
  --o-denoising-stats stats-dada2.qza

qiime metadata tabulate \
  --m-input-file stats-dada2.qza \
  --o-visualization stats-dada2.qzv
```
## Importing validated metadata
```
qiime metadata tabulate \
  --m-input-file v_metadata.tsv \
  --o-visualization tabulated-sample-metadata.qzv
```
## Generating Feature tables 
These show how many sequences are associated with each sample and with each feature, histograms of those distributions and some related summary statistics
```
qiime feature-table summarize \
  --i-table table-dada2.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file v_metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-dada2.qza \
  --o-visualization rep-seqs.qzv
```
## Generating phylogenetic trees
```
# Perform a multiple sequence alignment using the mafft program
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-dada2.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
```
## Rarefaction plotting
```
# Checking whether more sequencing depth could give different results
qiime diversity alpha-rarefaction \
  --i-table table-dada2.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 20000 \
  --m-metadata-file v_metadata.tsv \
  --o-visualization alpha-rarefaction.qzv
```
## Diversity analysis
Checking the diversity of each sample (how different are the individual samples and groups of samples)
### Alpha diversity
Shannon’s diversity index (a quantitative measure of community richness)
Observed Features (a qualitative measure of community richness)
Faith’s Phylogenetic Diversity (a qualitative measure of community richness that incorporates phylogenetic relationships between the features)
Evenness (or Pielou’s Evenness; a measure of community evenness)
### Beta diversity
Jaccard distance (a qualitative measure of community dissimilarity)
Bray-Curtis distance (a quantitative measure of community dissimilarity)
unweighted UniFrac distance (a qualitative measure of community dissimilarity that incorporates phylogenetic relationships between the features)
weighted UniFrac distance (a quantitative measure of community dissimilarity that incorporates phylogenetic relationships between the features)
```
# Options: --p-sampling-depth Choose a value that is as high as possible (so you retain more sequences per sample) while excluding as few samples as possible
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table-dada2.qza \
  --p-sampling-depth 20000 \
  --m-metadata-file v_metadata.tsv \
  --output-dir core-metrics-results
```
## Taxonomical analysis
Download the sequences from training from databases eg silva, greengenes, rmd and others
### Greengenes classification
Downloading the classifier
```
wget -C	` https://data.qiime2.org/classifiers/greengenes/gg_2022_10_backbone_full_length.nb.qza
```
## Runing the classification
```
qiime feature-classifier classify-sklearn \
  --i-classifier gg_2022_10_backbone_full_length.nb.qza \
  --i-reads rep-seqs-dada2.qza \
  --o-classification gg-taxonomy.qza

qiime metadata tabulate \
  --m-input-file gg-taxonomy.qza \
  --o-visualization gg-taxonomy.qzv
```
```
END="$( date +%s )"
echo "==========================This pipeline took: $[ $START - $END ] seconds==========================="
```
