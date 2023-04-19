###################################
#### preprocess the count data ####
###################################

rm(list=ls())
library(readr)
library(qiime2R)
library(phyloseq)
library(dplyr)

# set working directory to be where the current script is located
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

setwd('~/data/')

## read in the data
taxon_all <- read_qza('taxonomy.qza')$data # taxanomy of 3835 ASV from larger table
dim(taxon_all)
otu_tab <- read_qza('table.qza')$data # smaller ASV table 
dim(otu_tab)
meta_tab <- read_tsv('metadata-matched.tsv') # metadata for the smaller ASV table
dim(meta_tab)
colnames(meta_tab)[1] <- 'SampleID'

## retrieve taxonomy info for the smaller ASV table
taxon_tab <- as.matrix(parse_taxonomy(taxon_all))
ind <- match(rownames(otu_tab), rownames(taxon_tab))
taxon_tab <- taxon_tab[ind,]
table(rownames(otu_tab)==rownames(taxon_tab))

## format into phyloseq object
meta_tab <- as.data.frame(meta_tab)
rownames(meta_tab) <- meta_tab[,1]
table(colnames(otu_tab)==rownames(meta_tab))
ecam <- phyloseq(otu_table(otu_tab, taxa_are_rows = T), tax_table(taxon_tab), sample_data(meta_tab))
ecam_genus0 <- tax_glom(ecam, taxrank='Genus') # aggregate at genus level

## data cleaning
ecam_clean <- ecam_genus0
sample_variables(ecam_clean)
# remove samples from mother
ecam_clean <- subset_samples(ecam_clean, mom_child=='C')
# remove month 15 and 19
ecam_clean <- subset_samples(ecam_clean, !month %in% c(15,19))
# remove genus present in <5 samples
ecam_clean <- filter_taxa(ecam_clean, function(x) sum(x!=0) >= 5, TRUE)
# remove samples with <2000 reads
ecam_clean <- prune_samples(sample_sums(ecam_clean)>=2000, ecam_clean)
# remove subject with only one sample
tb <- with(sample_data(ecam_clean), table(month, studyid))
colnames(tb)[colSums(tb)<=1]
ecam_clean <- subset_samples(ecam_clean, studyid!=23)
tb <- with(sample_data(ecam_clean), table(month, studyid))
tb

ecam_clean
count_tab <- t(as.matrix(otu_table(ecam_clean)))
meta_tab <- as.matrix(sample_data(ecam_clean))
taxon_tab <- as.matrix(tax_table(ecam_clean))
write.csv(count_tab, file='genus_count_cleaned_q2.csv')
write.csv(meta_tab, file='genus_metadata_cleaned_q2.csv')
write.csv(taxon_tab, file='genus_taxonomy_cleaned_q2.csv')


ecam_clean <- ecam
sample_variables(ecam_clean)
# remove samples from mother
ecam_clean <- subset_samples(ecam_clean, mom_child=='C')
# remove month 15 and 19
ecam_clean <- subset_samples(ecam_clean, !month %in% c(15,19))
# remove genus present in <5 samples
ecam_clean <- filter_taxa(ecam_clean, function(x) sum(x!=0) >= 5, TRUE)
# remove samples with <2000 reads
ecam_clean <- prune_samples(sample_sums(ecam_clean)>=2000, ecam_clean)
# remove subject with only one sample
tb <- with(sample_data(ecam_clean), table(month, studyid))
colnames(tb)[colSums(tb)<=1]
ecam_clean <- subset_samples(ecam_clean, studyid!=23)
tb <- with(sample_data(ecam_clean), table(month, studyid))
tb

ecam_clean
count_tab <- t(as.matrix(otu_table(ecam_clean)))
meta_tab <- as.matrix(sample_data(ecam_clean))
taxon_tab <- as.matrix(tax_table(ecam_clean))
write.csv(count_tab, file='otu_count_cleaned_q2.csv')
write.csv(meta_tab, file='otu_metadata_cleaned_q2.csv')
write.csv(taxon_tab, file='otu_taxonomy_cleaned_q2.csv')


