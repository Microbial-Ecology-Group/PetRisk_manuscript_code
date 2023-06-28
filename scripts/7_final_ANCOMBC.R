library(ANCOMBC)
library(dplyr)
library(caret)
library(DT)
library(tidyr)
library(RColorBrewer)

# Rename the Group level taxonomy to facilitate aggregation with ANCOMBC
#colnames(phyloseq::tax_table(noSNP))[4] <- "Species"
#colnames(phyloseq::tax_table(noSNP))[2] <- "Class"
sample_data_df <- sample_data(noSNP)
lapply(sample_data_df, function(x) if (is.factor(x)) length(unique(x)))

sample_data(noSNP)$household_type <- as.factor(sample_data(noSNP)$household_type)
sample_data(noSNP)$location <- as.factor(sample_data(noSNP)$location)
sample_data(noSNP)$Group <- factor(sample_data(noSNP)$Group, levels = c("dog_healthy","dog_infected","human_healthy","human_infected"))


# All dogs
Dogs_noSNP <- subset_samples(noSNP, organism=="dog")
Dogs_group_noSNP <- tax_glom(Dogs_noSNP, taxrank = "group")
Dogs_group_noSNP <- prune_taxa(taxa_sums(Dogs_group_noSNP)>7, Dogs_group_noSNP)



# Just portugal
Portugal_noSNP <- subset_samples(noSNP, location=="Portugal")
Portugal_group_noSNP <- tax_glom(Portugal_noSNP, taxrank = "group")
Portugal_group_noSNP <- prune_taxa(taxa_sums(Portugal_group_noSNP)>7, Portugal_group_noSNP)


Humans_noSNP <- subset_samples(noSNP, organism=="human")
Humans_group_noSNP <- tax_glom(Humans_noSNP, taxrank = "group")
Humans_group_noSNP <- prune_taxa(taxa_sums(Humans_group_noSNP)>7, Humans_group_noSNP)

Portugal_Humans_noSNP <- subset_samples(noSNP, organism=="human" & location=="Portugal")
Portugal_Humans_group_noSNP <- tax_glom(Portugal_Humans_noSNP, taxrank = "group")
Portugal_Humans_group_noSNP <- prune_taxa(taxa_sums(Portugal_Humans_group_noSNP)>7, Portugal_Humans_group_noSNP)


Portugal_Dogs_noSNP <- subset_samples(noSNP, organism=="dog" & location=="Portugal")
Portugal_Dogs_group_noSNP <- tax_glom(Portugal_Dogs_noSNP, taxrank = "group")
Portugal_Dogs_group_noSNP <- prune_taxa(taxa_sums(Portugal_Dogs_group_noSNP)>7, Portugal_Dogs_group_noSNP)

Healthy_Dogs_noSNP <- subset_samples(noSNP, organism=="dog" & household_type == "healthy")
Healthy_Dogs_group_noSNP <- tax_glom(Healthy_Dogs_noSNP, taxrank = "group")
Healthy_Dogs_group_noSNP <- prune_taxa(taxa_sums(Healthy_Dogs_group_noSNP)>7, Healthy_Dogs_group_noSNP)



UK_group_noSNP <- tax_glom(UK_noSNP, taxrank = "group")
UK_group_noSNP <- prune_taxa(taxa_sums(UK_group_noSNP)>7, UK_group_noSNP)

microbiome_group.ps <- tax_glom(noSNP, taxrank = "group")



#### All dogs, by group ####
colnames(phyloseq::tax_table(Dogs_group_noSNP))[4] <- "Species"

output_group_dogs = ancombc2(data = Dogs_group_noSNP, tax_level = "Species",
                        fix_formula ="location + household_type", rand_formula = NULL,
                        p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE,
                        prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                        group = "household_type", struc_zero = TRUE, neg_lb = TRUE,
                        alpha = 0.05, n_cl = 2, verbose = TRUE,
                        global = TRUE, pairwise = TRUE, dunnet = FALSE, trend = FALSE,
                        iter_control = list(tol = 1e-2, max_iter = 20, 
                                            verbose = TRUE),
                        em_control = list(tol = 1e-5, max_iter = 20),
                        lme_control = lme4::lmerControl(),
                        mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                        trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                    nrow = 2, 
                                                                    byrow = TRUE),
                                                             matrix(c(-1, 0, 1, -1),
                                                                    nrow = 2, 
                                                                    byrow = TRUE)),
                                             node = list(2, 2),
                                             solver = "ECOS",
                                             B = 100))




## Structural zeros ####
tab_zero = output_group_dogs$zero_ind
tab_zero %>%
  datatable(caption = "The detection of structural zeros")

## Sensitivity scores ####
tab_sens = output_group_dogs$pseudo_sens_tab
tab_sens %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(tab_sens)[-1], digits = 2)



### ANCOM-BC2 primary analysis ####
#Result from the ANCOM-BC2 methodology to determine taxa that are differentially abundant
#according to the covariate of interest. Results contain: 1) log fold changes, 2) standard errors, 
#3) test statistics, 4) p-values, 5) adjusted p-values, 6) indicators of whether the taxon is
#differentially abundant (TRUE) or not (FALSE).

res_prim = output_group_dogs$res

## Results for household_type order  #####
# Original results
df_household = res_prim %>%
  dplyr::select(taxon, contains("household")) 


df_fig_household = df_household %>%
  filter(diff_household_typeinfected == 1) %>%
  mutate(lfc_household_type = ifelse(diff_household_typeinfected == 1, 
                                     lfc_household_typeinfected, 0)) %>%
  transmute(taxon, `Infected vs Healthy` = round(lfc_household_type, 2)) %>%
  pivot_longer(cols = `Infected vs Healthy`, names_to = "group", values_to = "value") %>%
  arrange(-value) %>%
  mutate(qvalue = df_household$q_household_typeinfected[match(taxon, df_household$taxon)])



fig_household = df_fig_household %>%
  ggplot(aes(x = group, y = reorder(taxon, value), fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradientn(colors =  rev(brewer.pal(9, "Blues")[c(1:8)]),
                       values = rescale(-c(up, lo), to = c(0, 1)),
                       na.value = "white", 
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "LogFC in dogs compared by household type") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_household


##
### Structural zeros ####
##
tab_zero = output_group$zero_ind
tab_zero %>%
  datatable(caption = "The detection of structural zeros")



# join the two dataframes by the taxon column
joined_df <- left_join(df_fig_household, tab_zero, by = "taxon")

# check if any significantly different taxons have structural zeros
any(!is.na(joined_df$`structural_zero (household_type = healthy)`) | 
      !is.na(joined_df$`structural_zero (household_type = infected)`))



# Get taxons with structural zeros for healthy household type
structural_zero_healthy <- tab_zero[tab_zero$`structural_zero (household_type = healthy)` == TRUE,]$taxon

# Get taxons with structural zeros for infected household type
structural_zero_infected <- tab_zero[tab_zero$`structural_zero (household_type = infected)` == TRUE,]$taxon

# Check if significantly different taxons match the structural zeros
match_healthy <- df_fig_household$taxon %in% structural_zero_healthy
match_infected <- df_fig_household$taxon %in% structural_zero_infected

# Print taxons that have structural zeros in the ANCOMBC2 model
cat("Taxons with structural zeros for healthy household type:\n")
print(structural_zero_healthy)
cat("\n")

cat("Taxons with structural zeros for infected household type:\n")
print(structural_zero_infected)
cat("\n")

# Print whether the significantly different taxons match the structural zeros
cat("Do significantly different taxons match the structural zeros for healthy household type?\n")
print(all(match_healthy == structural_zero_healthy))

cat("Do significantly different taxons match the structural zeros for infected household type?\n")
print(all(match_infected == structural_zero_infected))





##
### Sensitivity scores ####
##

tab_sens = output_group_dogs$pseudo_sens_tab
tab_sens %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(tab_sens)[-1], digits = 2)


sens_type = tab_sens %>%
  transmute(taxon, sens_type = household_typeinfected) %>%
  left_join(df_household_type %>%
              transmute(taxon, diff_type = diff_household_typeinfected), 
            by = "taxon") %>%
  mutate(group = "Infected vs. Healthy")

sens_type$diff_type = recode(sens_type$diff_type * 1, 
                             `1` = "Significant",
                             `0` = "Nonsignificant")

fig_sens_type = sens_type %>%
  ggplot(aes(x = taxon, y = sens_type, color = diff_type)) +
  geom_point() +
  scale_color_brewer(palette = "Dark2", name = NULL) +
  facet_grid(rows = vars(group), scales = "free") +
  labs(x = NULL, y = "Sensitivity Score") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5))
fig_sens_type






##
####
##### Just Portugal dogs #####
####
##



colnames(phyloseq::tax_table(Portugal_Dogs_group_noSNP))[4] <- "Species"

# Portugal_noSNP <- subset_samples(Dogs_noSNP, location=="Portugal")
# Portugal_noSNP <- prune_taxa(taxa_sums(Portugal_noSNP)>0, Portugal_noSNP)


output_group_portugal = ancombc2(data = Portugal_Dogs_group_noSNP , tax_level = "Species",
                        fix_formula ="household_type", rand_formula = NULL,
                        p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE,
                        prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                        group = "household_type", struc_zero = TRUE, neg_lb = TRUE,
                        alpha = 0.05, n_cl = 2, verbose = TRUE,
                        global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                        iter_control = list(tol = 1e-2, max_iter = 20, 
                                            verbose = TRUE),
                        em_control = list(tol = 1e-5, max_iter = 20),
                        lme_control = lme4::lmerControl(),
                        mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                        trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                    nrow = 2, 
                                                                    byrow = TRUE),
                                                             matrix(c(-1, 0, 1, -1),
                                                                    nrow = 2, 
                                                                    byrow = TRUE)),
                                             node = list(2, 2),
                                             solver = "ECOS",
                                             B = 100))


### ANCOM-BC2 primary analysis ####
#Result from the ANCOM-BC2 methodology to determine taxa that are differentially abundant
#according to the covariate of interest. Results contain: 1) log fold changes, 2) standard errors, 
#3) test statistics, 4) p-values, 5) adjusted p-values, 6) indicators of whether the taxon is
#differentially abundant (TRUE) or not (FALSE).

res_prim = output_group_portugal$res

## Results for household_type order  #####
# Original results
df_household = res_prim %>%
  dplyr::select(taxon, contains("household")) 


# df_fig_household = df_household %>%
#   filter(diff_household_typeinfected == 1 ) %>%
#   mutate(lfc_household_type = ifelse(diff_household_typeinfected == 1, 
#                                       lfc_household_typeinfected, 0)) %>%
#   transmute(taxon, 
#             `Infected vs Healthy` = round(lfc_household_type, 2)) %>%
#   pivot_longer(cols = `Infected vs Healthy`, 
#                names_to = "group", values_to = "value") %>%
#   arrange(-value)


df_fig_household = df_household %>%
  filter(diff_household_typeinfected == 1) %>%
  mutate(lfc_household_type = ifelse(diff_household_typeinfected == 1, 
                                     lfc_household_typeinfected, 0)) %>%
  transmute(taxon, `Infected vs Healthy` = round(lfc_household_type, 2)) %>%
  pivot_longer(cols = `Infected vs Healthy`, names_to = "group", values_to = "value") %>%
  arrange(-value) %>%
  mutate(qvalue = df_household$q_household_typeinfected[match(taxon, df_household$taxon)])





fig_household = df_fig_household %>%
  ggplot(aes(x = group, y = reorder(taxon, value), fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradientn(colors =  rev(brewer.pal(9, "Blues")[c(1:8)]),
                       values = rescale(-c(up, lo), to = c(0, 1)),
                       na.value = "white", 
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "LogFC in Portugal dogs compared by household type") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_household


##
### Structural zeros ####
##
tab_zero = output_group_portugal$zero_ind
tab_zero %>%
  datatable(caption = "The detection of structural zeros")



# join the two dataframes by the taxon column
joined_df <- left_join(df_fig_household, tab_zero, by = "taxon")

# check if any significantly different taxons have structural zeros
any(!is.na(joined_df$`structural_zero (household_type = healthy)`) | 
      !is.na(joined_df$`structural_zero (household_type = infected)`))



# Get taxons with structural zeros for healthy household type
structural_zero_healthy <- tab_zero[tab_zero$`structural_zero (household_type = healthy)` == TRUE,]$taxon

# Get taxons with structural zeros for infected household type
structural_zero_infected <- tab_zero[tab_zero$`structural_zero (household_type = infected)` == TRUE,]$taxon

# Check if significantly different taxons match the structural zeros
match_healthy <- df_fig_household$taxon %in% structural_zero_healthy
match_infected <- df_fig_household$taxon %in% structural_zero_infected

# Print taxons that have structural zeros in the ANCOMBC2 model
cat("Taxons with structural zeros for healthy household type:\n")
print(structural_zero_healthy)
cat("\n")

cat("Taxons with structural zeros for infected household type:\n")
print(structural_zero_infected)
cat("\n")

# Print whether the significantly different taxons match the structural zeros
cat("Do significantly different taxons match the structural zeros for healthy household type?\n")
print(all(match_healthy == structural_zero_healthy))

cat("Do significantly different taxons match the structural zeros for infected household type?\n")
print(all(match_infected == structural_zero_infected))





##
### Sensitivity scores ####
##

tab_sens = output_group_portugal$pseudo_sens_tab
tab_sens %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(tab_sens)[-1], digits = 2)


sens_type = tab_sens %>%
  transmute(taxon, sens_type = household_typeinfected) %>%
  left_join(df_household_type %>%
              transmute(taxon, diff_type = diff_household_typeinfected), 
            by = "taxon") %>%
  mutate(group = "Infected vs. Healthy")

sens_type$diff_type = recode(sens_type$diff_type * 1, 
                             `1` = "Significant",
                             `0` = "Nonsignificant")

fig_sens_type = sens_type %>%
  ggplot(aes(x = taxon, y = sens_type, color = diff_type)) +
  geom_point() +
  scale_color_brewer(palette = "Dark2", name = NULL) +
  facet_grid(rows = vars(group), scales = "free") +
  labs(x = NULL, y = "Sensitivity Score") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5))
fig_sens_type




##
####
##### Just humans #####
####
##



colnames(phyloseq::tax_table(Humans_group_noSNP))[4] <- "Species"

output_group_humans = ancombc2(data = Humans_group_noSNP , tax_level = "Species",
                        fix_formula ="household_type + location", rand_formula = NULL,
                        p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE,
                        prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                        group = "household_type", struc_zero = TRUE, neg_lb = TRUE,
                        alpha = 0.05, n_cl = 2, verbose = TRUE,
                        global = TRUE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                        iter_control = list(tol = 1e-2, max_iter = 20, 
                                            verbose = TRUE),
                        em_control = list(tol = 1e-5, max_iter = 20),
                        lme_control = lme4::lmerControl(),
                        mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                        trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                    nrow = 2, 
                                                                    byrow = TRUE),
                                                             matrix(c(-1, 0, 1, -1),
                                                                    nrow = 2, 
                                                                    byrow = TRUE)),
                                             node = list(2, 2),
                                             solver = "ECOS",
                                             B = 100))





### ANCOM-BC2 primary analysis ####
#Result from the ANCOM-BC2 methodology to determine taxa that are differentially abundant
#according to the covariate of interest. Results contain: 1) log fold changes, 2) standard errors, 
#3) test statistics, 4) p-values, 5) adjusted p-values, 6) indicators of whether the taxon is
#differentially abundant (TRUE) or not (FALSE).

res_prim = output_group_humans$res

## Results for household_type order  #####
# Original results
df_household = res_prim %>%
  dplyr::select(taxon, contains("household")) 


# df_fig_household = df_household %>%
#   filter(diff_household_typeinfected == 1 ) %>%
#   mutate(lfc_household_type = ifelse(diff_household_typeinfected == 1, 
#                                       lfc_household_typeinfected, 0)) %>%
#   transmute(taxon, 
#             `Infected vs Healthy` = round(lfc_household_type, 2)) %>%
#   pivot_longer(cols = `Infected vs Healthy`, 
#                names_to = "group", values_to = "value") %>%
#   arrange(-value)


df_fig_household = df_household %>%
  filter(diff_household_typeinfected == 1) %>%
  mutate(lfc_household_type = ifelse(diff_household_typeinfected == 1, 
                                     lfc_household_typeinfected, 0)) %>%
  transmute(taxon, `Infected vs Healthy` = round(lfc_household_type, 2)) %>%
  pivot_longer(cols = `Infected vs Healthy`, names_to = "group", values_to = "value") %>%
  arrange(-value) %>%
  mutate(qvalue = df_household$q_household_typeinfected[match(taxon, df_household$taxon)])


fig_household = df_fig_household %>%
  ggplot(aes(x = group, y = reorder(taxon, value), fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradientn(colors =  rev(brewer.pal(9, "Blues")[c(1:8)]),
                       values = rescale(-c(up, lo), to = c(0, 1)),
                       na.value = "white", 
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "LogFC in humans compared by household type") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_household


##
### Structural zeros ####
##
tab_zero = output_group_humans$zero_ind
tab_zero %>%
  datatable(caption = "The detection of structural zeros")



# join the two dataframes by the taxon column
joined_df <- left_join(df_fig_household, tab_zero, by = "taxon")

# check if any significantly different taxons have structural zeros
any(!is.na(joined_df$`structural_zero (household_type = healthy)`) | 
      !is.na(joined_df$`structural_zero (household_type = infected)`))



# Get taxons with structural zeros for healthy household type
structural_zero_healthy <- tab_zero[tab_zero$`structural_zero (household_type = healthy)` == TRUE,]$taxon

# Get taxons with structural zeros for infected household type
structural_zero_infected <- tab_zero[tab_zero$`structural_zero (household_type = infected)` == TRUE,]$taxon

# Check if significantly different taxons match the structural zeros
match_healthy <- df_fig_household$taxon %in% structural_zero_healthy
match_infected <- df_fig_household$taxon %in% structural_zero_infected

# Print taxons that have structural zeros in the ANCOMBC2 model
cat("Taxons with structural zeros for healthy household type:\n")
print(structural_zero_healthy)
cat("\n")

cat("Taxons with structural zeros for infected household type:\n")
print(structural_zero_infected)
cat("\n")

# Print whether the significantly different taxons match the structural zeros
cat("Do significantly different taxons match the structural zeros for healthy household type?\n")
print(all(match_healthy == structural_zero_healthy))

cat("Do significantly different taxons match the structural zeros for infected household type?\n")
print(all(match_infected == structural_zero_infected))





##
### Sensitivity scores ####
##

tab_sens = output_group_humans$pseudo_sens_tab
tab_sens %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(tab_sens)[-1], digits = 2)


sens_type = tab_sens %>%
  transmute(taxon, sens_type = household_typeinfected) %>%
  left_join(df_household_type %>%
              transmute(taxon, diff_type = diff_household_typeinfected), 
            by = "taxon") %>%
  mutate(group = "Infected vs. Healthy")

sens_type$diff_type = recode(sens_type$diff_type * 1, 
                             `1` = "Significant",
                             `0` = "Nonsignificant")

fig_sens_type = sens_type %>%
  ggplot(aes(x = taxon, y = sens_type, color = diff_type)) +
  geom_point() +
  scale_color_brewer(palette = "Dark2", name = NULL) +
  facet_grid(rows = vars(group), scales = "free") +
  labs(x = NULL, y = "Sensitivity Score") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5))
fig_sens_type


##
####
##### Just Portugal humans #####
####
##



colnames(phyloseq::tax_table(Portugal_Humans_group_noSNP))[4] <- "Species"


output_group_portugal_human = ancombc2(data = Portugal_Humans_group_noSNP , tax_level = "Species",
                        fix_formula ="household_type", rand_formula = NULL,
                        p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE,
                        prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                        group = "household_type", struc_zero = TRUE, neg_lb = TRUE,
                        alpha = 0.05, n_cl = 2, verbose = TRUE,
                        global = TRUE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                        iter_control = list(tol = 1e-2, max_iter = 20, 
                                            verbose = TRUE),
                        em_control = list(tol = 1e-5, max_iter = 20),
                        lme_control = lme4::lmerControl(),
                        mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                        trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                    nrow = 2, 
                                                                    byrow = TRUE),
                                                             matrix(c(-1, 0, 1, -1),
                                                                    nrow = 2, 
                                                                    byrow = TRUE)),
                                             node = list(2, 2),
                                             solver = "ECOS",
                                             B = 100))





### ANCOM-BC2 primary analysis ####
#Result from the ANCOM-BC2 methodology to determine taxa that are differentially abundant
#according to the covariate of interest. Results contain: 1) log fold changes, 2) standard errors, 
#3) test statistics, 4) p-values, 5) adjusted p-values, 6) indicators of whether the taxon is
#differentially abundant (TRUE) or not (FALSE).

res_prim = output_group_portugal_human$res

## Results for household_type order  #####
# Original results
df_household = res_prim %>%
  dplyr::select(taxon, contains("household")) 


# df_fig_household = df_household %>%
#   filter(diff_household_typeinfected == 1 ) %>%
#   mutate(lfc_household_type = ifelse(diff_household_typeinfected == 1, 
#                                       lfc_household_typeinfected, 0)) %>%
#   transmute(taxon, 
#             `Infected vs Healthy` = round(lfc_household_type, 2)) %>%
#   pivot_longer(cols = `Infected vs Healthy`, 
#                names_to = "group", values_to = "value") %>%
#   arrange(-value)


df_fig_household = df_household %>%
  filter(diff_household_typeinfected == 1) %>%
  mutate(lfc_household_type = ifelse(diff_household_typeinfected == 1, 
                                     lfc_household_typeinfected, 0)) %>%
  transmute(taxon, `Infected vs Healthy` = round(lfc_household_type, 2)) %>%
  pivot_longer(cols = `Infected vs Healthy`, names_to = "group", values_to = "value") %>%
  arrange(-value) %>%
  mutate(qvalue = df_household$q_household_typeinfected[match(taxon, df_household$taxon)])


fig_household = df_fig_household %>%
  ggplot(aes(x = group, y = reorder(taxon, value), fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradientn(colors =  rev(brewer.pal(9, "Blues")[c(1:8)]),
                       values = rescale(-c(up, lo), to = c(0, 1)),
                       na.value = "white", 
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "LogFC in humans compared by household type") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_household


##
### Structural zeros ####
##
tab_zero = output_group_portugal_human$zero_ind
tab_zero %>%
  datatable(caption = "The detection of structural zeros")



# join the two dataframes by the taxon column
joined_df <- left_join(df_fig_household, tab_zero, by = "taxon")

# check if any significantly different taxons have structural zeros
any(!is.na(joined_df$`structural_zero (household_type = healthy)`) | 
      !is.na(joined_df$`structural_zero (household_type = infected)`))



# Get taxons with structural zeros for healthy household type
structural_zero_healthy <- tab_zero[tab_zero$`structural_zero (household_type = healthy)` == TRUE,]$taxon

# Get taxons with structural zeros for infected household type
structural_zero_infected <- tab_zero[tab_zero$`structural_zero (household_type = infected)` == TRUE,]$taxon

# Check if significantly different taxons match the structural zeros
match_healthy <- df_fig_household$taxon %in% structural_zero_healthy
match_infected <- df_fig_household$taxon %in% structural_zero_infected

# Print taxons that have structural zeros in the ANCOMBC2 model
cat("Taxons with structural zeros for healthy household type:\n")
print(structural_zero_healthy)
cat("\n")

cat("Taxons with structural zeros for infected household type:\n")
print(structural_zero_infected)
cat("\n")

# Print whether the significantly different taxons match the structural zeros
cat("Do significantly different taxons match the structural zeros for healthy household type?\n")
print(all(match_healthy == structural_zero_healthy))

cat("Do significantly different taxons match the structural zeros for infected household type?\n")
print(all(match_infected == structural_zero_infected))





##
### Sensitivity scores ####
##

tab_sens = output_group_portugal_human$pseudo_sens_tab
tab_sens %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(tab_sens)[-1], digits = 2)


sens_type = tab_sens %>%
  transmute(taxon, sens_type = household_typeinfected) %>%
  left_join(df_household_type %>%
              transmute(taxon, diff_type = diff_household_typeinfected), 
            by = "taxon") %>%
  mutate(group = "Infected vs. Healthy")

sens_type$diff_type = recode(sens_type$diff_type * 1, 
                             `1` = "Significant",
                             `0` = "Nonsignificant")

fig_sens_type = sens_type %>%
  ggplot(aes(x = taxon, y = sens_type, color = diff_type)) +
  geom_point() +
  scale_color_brewer(palette = "Dark2", name = NULL) +
  facet_grid(rows = vars(group), scales = "free") +
  labs(x = NULL, y = "Sensitivity Score") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5))
fig_sens_type

##
####
##### Just healthy dogs #####
####
##



colnames(phyloseq::tax_table(Healthy_Dogs_group_noSNP))[4] <- "Species"


output_group_healthy_dogs = ancombc2(data = Healthy_Dogs_group_noSNP , tax_level = "Species",
                                       fix_formula ="location", rand_formula = NULL,
                                       p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE,
                                       prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                       group = "location", struc_zero = TRUE, neg_lb = TRUE,
                                       alpha = 0.05, n_cl = 2, verbose = TRUE,
                                       global = TRUE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                                       iter_control = list(tol = 1e-2, max_iter = 20, 
                                                           verbose = TRUE),
                                       em_control = list(tol = 1e-5, max_iter = 20),
                                       lme_control = lme4::lmerControl(),
                                       mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                       trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                                   nrow = 2, 
                                                                                   byrow = TRUE),
                                                                            matrix(c(-1, 0, 1, -1),
                                                                                   nrow = 2, 
                                                                                   byrow = TRUE)),
                                                            node = list(2, 2),
                                                            solver = "ECOS",
                                                            B = 100))





### ANCOM-BC2 primary analysis ####
#Result from the ANCOM-BC2 methodology to determine taxa that are differentially abundant
#according to the covariate of interest. Results contain: 1) log fold changes, 2) standard errors, 
#3) test statistics, 4) p-values, 5) adjusted p-values, 6) indicators of whether the taxon is
#differentially abundant (TRUE) or not (FALSE).

res_prim = output_group_healthy_dogs$res

## Results for household_type order  #####
# Original results
df_location = res_prim %>%
  dplyr::select(taxon, contains("location")) 

df_fig_location = df_location %>%
  filter(diff_locationUK == 1) %>%
  mutate(lfc_location = ifelse(diff_locationUK == 1, 
                                     lfc_locationUK, 0)) %>%
  transmute(taxon, `UK vs Portugal` = round(lfc_location, 2)) %>%
  pivot_longer(cols = `UK vs Portugal`, names_to = "group", values_to = "value") %>%
  arrange(-value) %>%
  mutate(qvalue = df_location$q_locationUK[match(taxon, df_location$taxon)])



fig_location = df_fig_location %>%
  ggplot(aes(x = group, y = reorder(taxon, value), fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradientn(colors =  brewer.pal(9, "Reds")[c(1:8)],
                       values = rescale(-c(up, lo), to = c(0, 1)),
                       na.value = "white", 
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "LogFC in healthy dogs compared by location") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_location


##
### Structural zeros ####
##
tab_zero = output_group_healthy_dogs$zero_ind
tab_zero %>%
  datatable(caption = "The detection of structural zeros")



# join the two dataframes by the taxon column
joined_df <- left_join(df_fig_location, tab_zero, by = "taxon")

# check if any significantly different taxons have structural zeros
any(!is.na(joined_df$`structural_zero (location = UK)`) | 
      !is.na(joined_df$`structural_zero (location = Portugal)`))



# Get taxons with structural zeros for healthy household type
structural_zero_UK <- tab_zero[tab_zero$`structural_zero (location = UK)` == TRUE,]$taxon

# Get taxons with structural zeros for infected household type
structural_zero_Portugal <- tab_zero[tab_zero$`structural_zero (location = Portugal)` == TRUE,]$taxon

# Check if significantly different taxons match the structural zeros
match_UK <- df_fig_location$taxon %in% structural_zero_UK
match_Portugal <- df_fig_location$taxon %in% structural_zero_Portugal



##
### Sensitivity scores ####
##

tab_sens = output_group_healthy_dogs$pseudo_sens_tab
tab_sens %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(tab_sens)[-1], digits = 2)


sens_type = tab_sens %>%
  transmute(taxon, sens_type = locationUK) %>%
  left_join(df_location %>%
              transmute(taxon, diff_type = diff_locationUK), 
            by = "taxon") %>%
  mutate(group = "UK vs Portugal")

sens_type$diff_type = recode(sens_type$diff_type * 1, 
                             `1` = "Significant",
                             `0` = "Nonsignificant")

fig_sens_type = sens_type %>%
  ggplot(aes(x = taxon, y = sens_type, color = diff_type)) +
  geom_point() +
  scale_color_brewer(palette = "Dark2", name = NULL) +
  facet_grid(rows = vars(group), scales = "free") +
  labs(x = NULL, y = "Sensitivity Score") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5))
fig_sens_type
