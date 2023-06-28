# Beta diversity
#############################################################################################
##############################         BETA DIVERSITY         #########################@#####
#############################################################################################
#############################################################################################

#### CSS transformation of counts
noSNP.css <- phyloseq_transform_css(noSNP, log = F)

#### create d.f. for beta-diversity metadata
noSNP.css.df <- as(sample_data(noSNP.css),"data.frame")

# ordinate it based on Bray-Curtis
noSNP.ord <- vegan::metaMDS(comm = t(otu_table(noSNP.css)), try = 20, trymax = 500, distance = "bray", autotransform = F)

# plot the ordination with 95% confidence ellipses coloured by groups
plot_ordination(noSNP.css, noSNP.ord, type = "samples", color = "Group_country") +
  theme_bw() +
  geom_point(size = 2.5) +
  stat_ellipse(geom = "polygon", aes(fill=Group_country), level = 0.95, lty =2, alpha= 0.3) +
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

## STATS BETWEEN ALL GROUPS
# create distance matrix we use to ordinate (Bray-Curtis)
noSNP.dist <- vegdist(t(otu_table(noSNP.css)), method = "bray")

## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
all.groups.permanova <- pairwise.adonis(noSNP.dist, noSNP.css.df$Group_country, perm = 9999, p.adjust.m = "BH")
all.groups.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
all.groups.disper <- betadisper(noSNP.dist, noSNP.css.df$Group_country)
all.groups.permdisp <- permutest(all.groups.disper, permutations = 9999, pairwise = T)
all.groups.permdisp # looks like a few are significant


#############################################################################################
##############################         BETA DIVERSITY - Group       #########################@#####
#############################################################################################
#############################################################################################

#### CSS transformation of counts
#noSNP.css <- phyloseq_transform_css(noSNP, log = F)
ps_group.css <- tax_glom(noSNP.css, taxrank = "group")


#### create d.f. for beta-diversity metadata
ps_group.css.df <- as(sample_data(ps_group.css),"data.frame")

# ordinate it based on Bray-Curtis
group.ord <- vegan::metaMDS(comm = t(otu_table(ps_group.css)), try = 20, trymax = 500, distance = "bray", autotransform = F)

# plot the ordination with 95% confidence ellipses coloured by groups
plot_ordination(ps_group.css, noSNP.ord, type = "samples", color = "Group_country") +
  theme_bw() +
  geom_point(size = 2.5) +
  stat_ellipse(geom = "polygon", aes(fill=Group_country), level = 0.95, lty =2, alpha= 0.3) +
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

## STATS BETWEEN ALL GROUPS
# create distance matrix we use to ordinate (Bray-Curtis)
noSNP.dist <- vegdist(t(otu_table(ps_group.css)), method = "bray")

## pairwise PERMANOVA with 9999 permutations and Benjamini-Hochberg correction
groups.permanova <- pairwise.adonis(noSNP.dist, ps_group.css.df$Group_country, perm = 9999, p.adjust.m = "BH")
groups.permanova # there are significant differences between some groups

## checking the dispersion of variance via PERMDISP to confirm they aren't significantly different
groups.disper <- betadisper(noSNP.dist, ps_group.css.df$Group_country)
groups.permdisp <- permutest(groups.disper, permutations = 9999, pairwise = T)
groups.permdisp # looks like a few are significant




