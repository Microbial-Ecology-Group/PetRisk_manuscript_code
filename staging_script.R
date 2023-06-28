# Staging script


# This script loads the data, removes genes requiring SNP confirmation, and removes the household "
# Output phyloseq object is "noSNP"
source("scripts/1_load_v3_data.R")

# Performs alpha diversity calculations, Wilcoxon testing, and boxplot creation
source("R_analysis/2_alpha_diversity.R")

# Performs beta diversity and ordination testing with PERMANOVA
source("R_analysis/3_beta_diversity.R")

# Relative abundance plots
source("R_analysis/4_relative_abundance_dendro.R")


# Output other useful files
write.csv(otu_table(noSNP),"../Writing/Results/SNPremoved_AMR_count_matrix.csv")


write.csv(otu_table(data_noSNP.css),"../Writing/Results/CSSnorm_SNPremoved_AMR_count_matrix.csv")


