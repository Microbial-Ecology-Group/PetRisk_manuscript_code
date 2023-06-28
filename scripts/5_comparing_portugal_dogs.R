# Subset samples


Portugal_dog_health_pre_infected <- subset_samples(data_noSNP, original_group == "Portugal_dog_healthy" | original_group == "Portugal_dog_infected_pre"  )
Portugal_dog_health_pre_infected <- prune_taxa(taxa_sums(Portugal_dog_health_pre_infected) > 0, Portugal_dog_health_pre_infected) 

###
####
##### Beta diversity #####
####
###

#### CSS transformation of counts
Portugal_dog_health_pre_infected.css <- phyloseq_transform_css(Portugal_dog_health_pre_infected, log = F)

#### create d.f. for beta-diversity metadata
Portugal_dog_health_pre_infected.css.df <- as(sample_data(Portugal_dog_health_pre_infected.css),"data.frame")

# ordinate it based on Bray-Curtis
Portugal_dog_health_pre_infected.css.ord <- vegan::metaMDS(comm = t(otu_table(Portugal_dog_health_pre_infected.css)), try = 20, trymax = 500, distance = "bray", autotransform = F)

# plot the ordination with 95% confidence ellipses coloured by groups
plot_comparing_portugal_dogs <- plot_ordination(Portugal_dog_health_pre_infected.css, Portugal_dog_health_pre_infected.css.ord, type = "samples", color = "original_group") +
  theme_bw() +
  geom_point(size = 2.5) +
  stat_ellipse(geom = "polygon", aes(fill=original_group), level = 0.95, lty =2, alpha= 0.3) +
  #scale_fill_manual(values = group_palette) +
  #scale_colour_manual(values = group_palette) +
  theme(plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        panel.border = element_rect(colour = "black", size =0.75),
        legend.key.size = unit(2, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 24, vjust = 1.75),
        axis.title.x = element_text(size = 24, vjust = -1.5))
plot_comparing_portugal_dogs

## STATS BETWEEN ALL GROUPS
# create distance matrix we use to ordinate (Bray-Curtis)
Portugal_dog_health_pre_infected.dist <- vegdist(t(otu_table(Portugal_dog_health_pre_infected.css)), method = "bray")

## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
Portugal_dog_health_pre_infected.permanova <- pairwise.adonis(Portugal_dog_health_pre_infected.dist, Portugal_dog_health_pre_infected.css.df$original_group, perm = 9999, p.adjust.m = "BH")
Portugal_dog_health_pre_infected.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
Portugal_dog_health_pre_infected.disper <- betadisper(Portugal_dog_health_pre_infected.dist, Portugal_dog_health_pre_infected.css.df$original_group)
Portugal_dog_health_pre_infected.permdisp <- permutest(Portugal_dog_health_pre_infected.disper, permutations = 9999, pairwise = T)
Portugal_dog_health_pre_infected.permdisp # 



####
####  Comparing humans 
####
Portugal_human_health_pre_infected <- subset_samples(data_noSNP,original_group == "Portugal_human_healthy" | original_group == "Portugal_human_infected_pre"  )
Portugal_human_health_pre_infected <- prune_taxa(taxa_sums(Portugal_human_health_pre_infected) > 0, Portugal_human_health_pre_infected) 

###
####
##### Beta diversity #####
####
###

#### CSS transformation of counts
Portugal_human_health_pre_infected.css <- phyloseq_transform_css(Portugal_human_health_pre_infected, log = F)

#### create d.f. for beta-diversity metadata
Portugal_human_health_pre_infected.css.df <- as(sample_data(Portugal_human_health_pre_infected.css),"data.frame")

# ordinate it based on Bray-Curtis
Portugal_human_health_pre_infected.css.ord <- vegan::metaMDS(comm = t(otu_table(Portugal_human_health_pre_infected.css)), try = 20, trymax = 500, distance = "bray", autotransform = F)

# plot the ordination with 95% confidence ellipses coloured by groups
plot_comparing_portugal_dogs <- plot_ordination(Portugal_human_health_pre_infected.css, Portugal_human_health_pre_infected.css.ord, type = "samples", color = "original_group") +
  theme_bw() +
  geom_point(size = 2.5) +
  stat_ellipse(geom = "polygon", aes(fill=original_group), level = 0.95, lty =2, alpha= 0.3) +
  #scale_fill_manual(values = group_palette) +
  #scale_colour_manual(values = group_palette) +
  theme(plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        panel.border = element_rect(colour = "black", size =0.75),
        legend.key.size = unit(2, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 24, vjust = 1.75),
        axis.title.x = element_text(size = 24, vjust = -1.5))
plot_comparing_portugal_dogs

## STATS BETWEEN ALL GROUPS
# create distance matrix we use to ordinate (Bray-Curtis)
Portugal_human_health_pre_infected.dist <- vegdist(t(otu_table(Portugal_human_health_pre_infected.css)), method = "bray")

## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
Portugal_human_health_pre_infected.permanova <- pairwise.adonis(Portugal_human_health_pre_infected.dist, Portugal_human_health_pre_infected.css.df$original_group, perm = 9999, p.adjust.m = "BH")
Portugal_human_health_pre_infected.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
Portugal_human_health_pre_infected.disper <- betadisper(Portugal_human_health_pre_infected.dist, Portugal_human_health_pre_infected.css.df$original_group)
Portugal_human_health_pre_infected.permdisp <- permutest(Portugal_human_health_pre_infected.disper, permutations = 9999, pairwise = T)
Portugal_human_health_pre_infected.permdisp # 


