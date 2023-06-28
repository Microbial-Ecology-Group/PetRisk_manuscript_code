# alpha diversity


#############################################################################################
##############################         ALPHA DIVERSITY         ##############################
#############################################################################################
#############################################################################################

#### calculate richness, Shannon, Simpson, and inverse Simpson
alpha_div <- estimate_richness(noSNP, measures = c("Observed","Shannon","Simpson","InvSimpson"))
alpha_div.df <- as(sample_data(noSNP), "data.frame") # creating a df of the metadata
alpha_div_meta <- cbind(alpha_div, alpha_div.df) # combinging the alpha div values and metadata

write.csv(alpha_div_meta,"figures/All_v3_diversity_indices.csv")

##$$ splitting up based on healthy/infected for within-group comparisons
#alpha_div_healthy <- alpha_div_meta[which(alpha_div_meta$household_type=="healthy"),]
#alpha_div_infected <- alpha_div_meta[which(alpha_div_meta$household_type=="infected"),]

##### boxplot plot for richness comparing across the groups
all_observed_by_group <- ggplot(alpha_div_meta, aes(x= Group_country, y= Observed, fill = Group_country, colour = Group_country)) +
  theme_bw() +
  labs(title="RICHNESS", y= "Observed AMR genes") +
  geom_boxplot(alpha = 0.65, size = 1) +
  #scale_fill_manual(values = group_palette) +
  #scale_colour_manual(values = group_palette) +
  theme(plot.margin = unit(c(0.5,1,1,1), "cm"),
        plot.title = element_text(size = 28),
        legend.key.size = unit(2, "lines"),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 28, vjust = 1.75),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
all_observed_by_group

#### pairwise Wilxocon rank-sum with Benjamini-Hochberg correction for mult. comps.
richness_groups.pw <- pairwise.wilcox.test(alpha_div_meta$Observed, alpha_div_meta$Group_country, p.adjust.method = "BH")
richness_groups.pw # p-values in the matrix






# Shannon
all_shannon_by_group <- ggplot(alpha_div_meta, aes(x= Group_country, y= Shannon, fill = Group_country, colour = Group_country)) +
  theme_bw() +
  labs(title="Shannon", y= "Shannon") +
  geom_boxplot(alpha = 0.65, size = 1) +
  #scale_fill_manual(values = group_palette) +
  #scale_colour_manual(values = group_palette) +
  theme(plot.margin = unit(c(0.5,1,1,1), "cm"),
        plot.title = element_text(size = 28),
        legend.key.size = unit(2, "lines"),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 28, vjust = 1.75),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
all_shannon_by_group


shannon_groups.pw <- pairwise.wilcox.test(alpha_div_meta$Shannon, alpha_div_meta$Group_country, p.adjust.method = "BH")
shannon_groups.pw # p-values in the matrix










