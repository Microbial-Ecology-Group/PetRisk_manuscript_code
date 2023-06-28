# Upset plots

#Taras Plots# Much simpler than Lee's but his still works

library(UpSetR)

library("MicrobiotaProcess")


group_noSNP.css <- tax_glom(noSNP.css, taxrank = "group")


upsetda_edit <- get_upset(group_noSNP.css, factorNames="Group_country") ## ASV
upset(upsetda_edit, sets=c("Portugal_human_healthy", "Portugal_human_infected", "Portugal_dog_healthy", "Portugal_dog_infected","UK_human_healthy","UK_dog_healthy"),

sets.bar.color = c("#d6e2e9","#e9d556","#ac1d1c", "blue","red","green"),text.scale = 2,

order.by = "freq", empty.intersections = "on")



# Upset plot for UK 

UK_group_noSNP.css <- subset_samples(group_noSNP.css, location=="UK")

upsetda_edit <- get_upset(UK_group_noSNP.css, factorNames="Group_country") ## ASV
upset(upsetda_edit, sets=c("UK_human_healthy","UK_dog_healthy"),
      
      sets.bar.color = c("green","red"),text.scale = 2,
      
      order.by = "freq", empty.intersections = "on")

##
## Upset plot for Portugal ####
##

Portugal_group_noSNP.css <- subset_samples(group_noSNP.css, location=="Portugal")

upsetda_edit <- get_upset(Portugal_group_noSNP.css, factorNames="Group_country") ## ASV
upset(upsetda_edit, sets=c("Portugal_human_healthy", "Portugal_human_infected", "Portugal_dog_healthy", "Portugal_dog_infected"),
      
      sets.bar.color = c("#d6e2e9", "blue","red","green"),text.scale = 2,
      
      order.by = "freq", empty.intersections = "on")
