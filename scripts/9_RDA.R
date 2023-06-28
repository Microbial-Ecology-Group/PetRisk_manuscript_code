
cleanset <- data.frame(sample_data(noSNP))
cleanset <- cleanset[, c("location","organism","household_type","household_number")]

############################################### RDA MODEL ###################################################

#### CSS transformation of counts
noSNP.css <- phyloseq_transform_css(noSNP, log = F)
group_noSNP.css <- tax_glom(noSNP.css, taxrank = "group")


norm.mt<-as(otu_table(group_noSNP.css),"matrix")
### This matrix has rows as otus and columns as samples (all our work before used columns as otus)
## we first remove all rows (otus) that have count zero by summing per row and keeping those rows with
## sum greater than zero
norm.rs = apply(norm.mt,1,sum)
norm.mt = norm.mt[norm.rs>0,]
## We then transform the normalized matrix so that rows are samples and columns are otus (for vegan to work correctly)
norm.mt = t(norm.mt)
#norm.mt<-norm.mt[cleanset_samples,]
#Using function decostand() in vegan to transform using hellinger standardizations
hell.norm.mt <- decostand(norm.mt,"hell")
### We define the full model that includes all possible variables in the metadata file
mod1 <- rda(hell.norm.mt ~ ., data =  cleanset)
anova(mod1, by = "term", perm = 1000)

### We define the null model that only includes an intercept (1)
mod0 <- rda(hell.norm.mt ~ 1, data =  cleanset)
### We use the step function to find the best model by removing one term at a time.
### We use permutation test (perMANOVA) to do so
mod <- step(mod0, scope = formula(mod1), test = "perm", direction="both")
anova(mod, by = "term", perm = 1000)
anova(mod,perm = 1000)

# View the anova results
mod

# Visualize
#### 2D RDA plot (this doesn't work for db-RDA, plots are meaningless there)
fig <-ordiplot(mod,type="n")
colvec <- c("black","grey51","blue","red","green") 
#pchvec <-c(10,2)
points(mod, "sites", pch=21, col=as.factor(scaleddf$var), bg=as.factor(scaleddf$var))
groupz <- sort(unique(scaleddf$var))
for(i in seq(groupz)) {ordispider(mod,  scaleddf$catvar,font=2, cex=1.5, col=i, show.groups=groupz[i])}
legv <- sort(unique(scaleddf$catvar))
#legend("right", legend = levels(scaled_metadata$Group),col=legv, title = "Sample groups",bty = "n", pt.bg = legv)
text(fig, "biplot", col="red", cex=0.9)
text(fig, "centroids", col=c("black","red","green"), cex=2, pos = 4)


# Apriori model
# We can define our own model here
mod_apriori <- rda(hell.norm.mt ~ contvar  + catvars , data =  scaleddf)
anova(mod_apriori,perm = 10000)
anova(mod_apriori, by = "term", perm = 10000)

fig <-ordiplot(mod_apriori,type="n")
fig <-ordiplot(mod_apriori,type="n")
colvec <- c("black","grey51","blue","red","green") 


