

###
####
##### Differential abundance #####
####
###



## using the DA.zig to run metagenomeSeq's fitZIG to look for DA AMR Classes

## agglomerate our CSS normalized count table at the class level
class.css <- tax_glom(noSNP.css, taxrank = "Class")

# only want to include classes making up greater than 0.1% mean relative abundance
filt_ra_class <- filter_taxa(rel_abund_class, function(x) mean(x) > 0.1, T) # leave us with 30 classes
taxa_to_keep <- taxa_names(tax_table(filt_ra_class))
taxa_to_keep

## remove low abundance features
filt_class.css <- prune_taxa(taxa_to_keep, class.css) 

#ZIG model for DA classes between organisms, with location and healthy/infected as covariables
DA_zig_organism_all <- DA.zig(filt_class.css, predictor = "organism", covars = c("location","household_type"), p.adj = "BH", allResults = F)

# view the results
DA_zig_organism_all

# write a table of the results
write.csv(DA_zig_organism_all, "ZIGfit_results_organism_all.csv")

### from the adj. p-values in that table there were quite a few classes deferentially abundant:
# human > dog: Mupirocin, Glycopeptides, Nucleosides, Penicol, MLS, Trimethoprim, Aminoglycosides
# dog > human: Chromium resistance, Bacitracin, Iron resistance

# lets just plot one to see what it looks like with RA (it's not an apple to apples comparison since it won't have covariables)
aminoglycosides_ra_class <- subset_taxa(rel_abund_class, Class=="Aminoglycosides") %>%
  psmelt()
ggplot(aminoglycosides_ra_class, aes(x= organism, y= Abundance, fill = organism)) +
  theme_bw() +
  facet_grid(~location) +
  geom_bar(stat = "summary", colour = "black") + 
  geom_errorbar(stat = "summary")
# they are generally more abundant in Portugal but the trend of being higher in humans vs dogs is present in both countries
# makes sense and agrees with the ZIG model finding

### second ZIG model looking within just infected households and comparing time points
infected_class.css <- subset_samples(class.css, household_type=="infected") # subset out just infected households
filt_infected_class.css <- prune_taxa(taxa_to_keep, infected_class.css) # keep taxa over 0.1% mean relative abundance only
sum(taxa_sums(filt_infected_class.css)==0) # check if there any taxa with 0 counts across all samples (nope, good to go)

#ZIG model comparing pre vs during and using organism as a covariable
DA_zig_infected_timepoint <- DA.zig(filt_infected_class.css, predictor = "time", covars = "organism", p.adj = "BH", allResults = F)

# view the results
DA_zig_infected_timepoint

# write the results
write.csv(DA_zig_infected_timepoint, "ZIGfit_results_infected_timepoint.csv")
# based on the table, no significantly deferentially abundant classes (biggest log fold-change was -0.241 in Glycopeptides)
# my guess is that there may be different trends within humans and dogs which cancels out any difference, 
#let's just take a look at Glycopeptides

glyco_ra_infected <- subset_taxa(infected_class.css, Class=="Glycopeptides") %>%
  psmelt()

# plot
ggplot(glyco_ra_infected, aes(x= time, y= Abundance, fill= time)) +
  theme_bw() +
  facet_grid(~organism) +
  geom_bar(stat = "summary", colour = "black") +
  geom_errorbar(stat = "summary")
# yep, goes down from pre-during in dogs but up in humans
# might want to split dogs from humans and try a new ZIG model, but all the tools you need for that are in script already!

## Good luck, feel free to ask me if you need help with anything
