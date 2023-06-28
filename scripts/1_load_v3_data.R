library(phyloseq);library(metagenomeSeq);library(dplyr);library(scales);
library(pairwiseAdonis); library(vegan); library(metagMisc); library(stringr)
library(ggplot2);  library(randomcoloR); library(cowplot)
library(pairwiseAdonis); library(picante); library(grid); library(gridExtra); 
library(ggalt); library(ggforce); library(concaveman); library(ggdendro)

# #### TE DATA ####
# read in AMR count matrix
amr_data <- read.table('bioinformatics/AMR_analytic_matrix.csv', header=T, row.names=1, sep=',', quote = "")
# convert this into an 'otu_table' format required for phyloseq
amr_data <- otu_table(amr_data, taxa_are_rows = T)

#read in gene annotations
annotations <- read.table('bioinformatics/megares_annotations_v3.00.csv', header=T, row.names=1, sep=",", quote = "")
#convert this into a 'taxonomy table' format required for phyloseq
annotations <- phyloseq::tax_table(as.matrix(annotations))

# read in TE metadata
amr_metadata <- read.table('bioinformatics/WeeseTE_metadata.txt', header=T, sep='\t', row.names = 1, quote = "")
# convert to 'sample_data' format for phyloseq
amr_metadata <- sample_data(amr_metadata)

# merge the annotations, the count matrix, and metadata into a phyloseq object
data <- merge_phyloseq(amr_data, annotations, amr_metadata)
#
# #### splitting of the genes requiring SNP confirmation
SNPconfirm <- subset_taxa(data, snp=="RequiresSNPConfirmation")
SNPconfirm # 284 of the 3054 genes require SNP confirmation
noSNP <- subset_taxa(data, snp!="RequiresSNPConfirmation")
noSNP # 2814 genes remain, these are what we will use from now on

# Remove household 118
noSNP <- subset_samples(noSNP, household_number != "P118") # removes 4 samples from the 118 household

#
# #### some QC checks
# any taxa with no reads?
sum(taxa_sums(noSNP)==0) # 0 taxa without any counts 
noSNP <- prune_taxa(taxa_sums(noSNP) > 0, noSNP) # redundant, but another check to remove taxa with 0 counts
# check the number of genes in the samples
sort(sample_sums(noSNP)) # 100,555 is lowest - 1,132,295 is highest; no problem there

