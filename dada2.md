---
title: "dada2"
author: "Syrus"
date: "2024-05-13"
output:
  pdf_document: default
  html_document: default
---
We used dada2 version 1.18 for this analysis. A total of 30 demultiplexed samples 

## Installing tools

```{r}
install.packages("BiocManager")
BiocManager::install("ShortRead", version = "3.12") 
BiocManager::install("additional_missing_package_1", version = "3.12") 
BiocManager::install("dada2")
```

## Setting working directory and importing files
```{r}
library(dada2); packageVersion("dada2")
setwd("C:/Users/endawula/Downloads/QIIME_ANALYSIS/QIIME_ANALYSIS/")
path <- "Microbiome_TB-Data/FASTQ_Files/"
list.files(path)
#list.files()

```

```{r}

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq.gz and 
# SAMPLENAME_R2.fastq.gz
fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XX.fastq.gz
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
View(fnFs)
```

## Inspect read quality profiles

```{r}
plotQualityProfile(fnFs, aggregate = T)
```


```{r}
plotQualityProfile(fnRs, aggregate = T)
```

In gray-scale is a heat map of the frequency of each quality score at each base position. The mean quality score at each position is shown by the green line, and the quartiles of the quality score distribution by the orange lines. The red line shows the scaled proportion of reads that extend to at least that position (this is more useful for other sequencing technologies, as Illumina reads are typically all the same length, hence the flat red line).

## Filter and trim
```{r}
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```
NB: Reads must still overlap after truncation in order to merge them later!

```{r}
FWD <- "GTGYCAGCMGCCGCGGTAA"
REV <- "GGACTACNVGGGTWTCTAAT"
trimLeft = c(FWD, REV)

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(300,260),
              maxN=0, maxEE = c(2,2), rm.phix=TRUE, compress=TRUE, 
              trimLeft = c(19, 20), multithread=F) # On Windows set multithread=FALSE
out
```
fnFs & fnRs - directory with forward and reverse reads respectively 
filtFs & filtRs - directories with filtered and trimmed sequences 
truncLen - Chops reads after a given base, reads shorter than this are discarded 
maxN - Removes sequences with Ns after truncation 
truncQ - Truncates reads at the first instance of a quality score less than integer
rm.phix - Discards reads that match against the phiX genome
minQ - Reads with quality below after truncation are discarded
maxEE - After truncation higher than expected error will be discarded
compress - When true output fastq files are gzipped
multithread - False for windows because doesnt support folking (mc.cores =1)

NB: Some of these parameters were left to default.

## Learn about error rates
The error rates for each possible transition (A→C, A→G, …) are shown. 
Points are the observed error rates for each consensus quality score. 
The black line shows the estimated error rates after convergence of the machine-learning algorithm. The red line shows the error rates expected under the nominal definition of the Q-score. Here the estimated error rates (black line) are a good fit to the observed rates (points), and the error rates drop with increased quality as expected. Everything looks reasonable and we proceed with confidence.
```{r}
errF <- learnErrors(filtFs, multithread=TRUE)

errR <- learnErrors(filtRs, multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
```
## Dereplication
Dereplication conbines all identical sequencing reads into "unique sequences" with a corresponding "abundance" equal to the number of reads with that unique sequence. This reduces computational time by eliminating redundant comparisons. It also retains a summary of the quality information associated with each unique sequence. The consensus quality profile of a unique is the average of the 
positional qualities from the dereplicated reads. These inform the error model of the subsequent sample inference step, significantly increasing DADA2's accuracy.
```{r}
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)

```

## Sample inference
```{r}
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
# Inspecting the returned dada-class object
dadaFs[[1]]
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
# Inspecting the returned dada-class object
dadaRs[[1]]

```

The DADA2 algorithm inferred 115 and 113 true sequence variants from the 3679 and 3994 unique sequences in the first forward and reverse sample respectively. 

## Merge paired reads
We now merge the forward and reverse reads together to obtain the full denoised sequences. Merging is performed by aligning the denoised forward reads with the reverse-complement of the corresponding denoised reverse reads, and then constructing the merged “contig” sequences. By default, merged sequences are only output if the forward and reverse reads overlap by at least 12 bases, and are identical to each other in the overlap region (but these conditions can be changed via function arguments).

```{r}
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
```
The mergers object is a list of data.frames from each sample. Each data.frame contains the merged $sequence, its $abundance, and the indices of the $forward and $reverse sequence variants that were merged. Paired reads that did not exactly overlap were removed by mergePairs, further reducing spurious output.

## Construct sequence table
We can now construct an amplicon sequence variant table (ASV) table, a higher-resolution version of the OTU table produced by traditional methods.
```{r}
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```
```{r}
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

```
The sequence table is a matrix with rows corresponding to (and named by) the samples, and columns corresponding to (and named by) the sequence variants. This table contains 4100 ASVs.

## Remove chimeras
The core dada method corrects substitution and indel errors, but chimeras remain. Fortunately, the accuracy of sequence variants after denoising makes identifying chimeric ASVs simpler than when dealing with fuzzy OTUs. Chimeric sequences are identified if they can be exactly reconstructed by combining a left-segment and a right-segment from two more abundant “parent” sequences.

```{r}
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
```

```{r}
sum(seqtab.nochim)/sum(seqtab)

```
The frequency of chimeric sequences varies substantially from dataset to dataset, and depends on factors including experimental procedures and sample complexity. Here chimeras make up about 25% of the merged sequence variants, but when we account for the abundances of those variants we see they account for only about 15% of the merged sequence reads.

## Track reads through the pipeline

```{r}
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), 
               sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: 
# e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", 
                     "merged", "nonchim")
rownames(track) <- sample.names
track

```

Looks good! We kept the majority of our raw reads, and there is no over-large drop associated with any single step. If a majority of reads failed to merge, you may need to revisit the truncLen parameter used in the filtering step and make sure that the truncated reads span your amplicon. If a majority of reads were removed as chimeric, you may need to revisit the removal of primers, as the ambiguous nucleotides in unremoved primers interfere with chimera identification

## Assign taxonomy

The DADA2 package provides a native implementation of the naive Bayesian classifier method for this purpose. The assignTaxonomy function takes as input a set of sequences to be classified and a training set of reference sequences 
with known taxonomy, and outputs taxonomic assignments with at least minBoot bootstrap confidence.

```{r}

taxa <- assignTaxonomy(seqtab.nochim, multithread=TRUE,  
                       "C:/Users/endawula/Downloads/GTDB_bac120_arc122_ssu_r202_fullTaxo.fa.gz")

```

Let’s inspect the taxonomic assignments.
```{r}
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```

Unsurprisingly, the Bacteroidetes are well represented among the most abundant taxa in these fecal samples. Few species assignments were made, both because it is often not possible to make unambiguous species assignments from subsegments of the 16S gene, and because there is surprisingly little coverage of the indigenous mouse gut microbiota in reference databases.

## Handoff to phyloseq

```{r}
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

```

We can construct a simple sample data.frame from the information encoded in the filenames. Usually this step would instead involve reading the sample data in from a file.

```{r}
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out

```

We now construct a phyloseq object directly from the dada2 outputs.

```{r}
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
```

It is more convenient to use short names for our ASVs (e.g. ASV21) rather than the full DNA sequence when working with some of the tables and visualizations from phyloseq, but we want to keep the full DNA sequences for other purposes like merging with other datasets or indexing into reference databases like the Earth Microbiome Project. For that reason we’ll store the DNA sequences of our ASVs in the refseq slot of the phyloseq object, and then rename our taxa to a short string. That way, the short new taxa names will appear in tables and plots, and we can still recover the DNA sequences corresponding to each ASV as needed with refseq(ps).

```{r}
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
```

We are now ready to use phyloseq!

### Visualise alpha-diversity

```{r}
plot_richness(ps, measures=c("Shannon", "Simpson"))

```

No obvious systematic difference in alpha-diversity between early and late samples.

```{r}
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

```

```{r}
plot_ordination(ps.prop, ord.nmds.bray, title="Bray NMDS")

```
Ordination picks out a clear separation between the early and late samples.

### Bar plot:
```{r}
top50 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:50]
ps.top50 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top50 <- prune_taxa(top50, ps.top50)
plot_bar(ps.top50, fill="Genus") + facet_wrap(~When, scales="free_x")

```

Nothing glaringly obvious jumps out from the taxonomic distribution of the top 50 sequences to explain the early-late differentiation.

These were minimal examples of what can be done with phyloseq, as our purpose here was just to show how the results of DADA2 can be easily imported into phyloseq and interrogated further. For examples of the many analyses possible with phyloseq, see the phyloseq web site!


    Here ends the DADA2 portion of the tutorial.