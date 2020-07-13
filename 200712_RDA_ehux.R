setwd("/Users/marianacamaradosreis/Documents/PhD/Chapter1/Final_results/analysis_R/200331_microbiome_analysis")

##libraries
library(vegan)
library(ggplot2)
library(lme4)
library(car)
library(ggpubr)
library (adespatial)
library(ggbiplot)

##Using table that came from pre-processing long_p
group <- read_tsv("200705_core_ehux_merged_newflag.txt")
print(group)
dim(group)
sum(group[,2:458])

#1846324
##1)Filtering out ASVs with less than 5% of the reads to plot in the heatmap##
otu <- column_to_rownames(group, var = "X1") %>%
  decostand(., method="hellinger") %>%   rownames_to_column(., "X1") %>%
  gather(var, value, -X1) %>% 
  filter(value >=0.2236) %>% pivot_wider(names_from = var, values_from = value) %>% 
  replace(., is.na(.), 0) %>% column_to_rownames(., var = "X1")

dim(otu)

sample_data <- read_tsv("200712_metadata_strains.txt") # each "c" means a column of "character"
print(sample_data)
dim(sample_data)

core_metadata <- otu %>%
  rownames_to_column(., "SampleID") %>%
  left_join(., sample_data, by = c("SampleID" = "sample")) %>%
  column_to_rownames(., "SampleID")

##2)Plotting the heatmap

plot <- pheatmap(core_metadata[,1:91], fontsize_col =7, fontsize_row = 8, color=colorRampPalette(rev(brewer.pal(n = 7, name = "Spectral")))(100), cluster_cols = FALSE, clustering_distance_rows = "euclidean", clustering_method="complete", cex=1, cutree_rows = 3)
     
##3) Extracting groups formed with clustering this I'll use to color the PCA

groups <- sort(cutree(plot$tree_row, k=3)) %>% as.data.frame(.) %>%
  rownames_to_column(., "strain")

write.table(groups, "200713_groups_5%.txt", sep="\t")

##4) PCA to show the main groups of microbiomes###
otu <- column_to_rownames(group, var = "X1") %>%
  decostand(., method="hellinger") 

core_metadata2 <- otu %>%
  rownames_to_column(., "SampleID") %>%
  left_join(., sample_data, by = c("SampleID" = "sample")) %>%
  column_to_rownames(., "SampleID")

dim(core_metadata2)

core.pca <- prcomp(core_metadata2[,1:457])

ggbiplot(core.pca)

ggbiplot(core.pca, var.axes=TRUE, groups=core_metadata2$groups,  circle=FALSE, ellipse=TRUE)


##5) Doing a heatmap with all the ASVs 
###
plot <- pheatmap(core_metadata2[,1:457], fontsize_col =3, fontsize_row = 8, color=colorRampPalette(rev(brewer.pal(n = 7, name = "Spectral")))(100), cluster_rows = TRUE, clustering_distance_rows = "euclidean", clustering_distance_cols = "correlation", clustering_method="ward.D", cex=1, cutree_rows = 3)



##6) Checking beta-dispersion of the groups for RDA##
###Checking betadispersion of the data in the groups####
dist_v = vegdist(core_metadata2[,1:457], method="euclidean")

disper_morphogroup=betadisper(dist_v, core_metadata2$morphogroup, type="centroid")
permutest(disper_morphogroup, permutations = 1000)
boxplot(disper_morphogroup)

disper_genus=betadisper(dist_v, core_metadata2$genus, type="centroid")
permutest(disper_genus)
boxplot(disper_genus)

disper_ocean=betadisper(dist_v, core_metadata2$ocean, type="centroid")
permutest(disper_ocean, control=permControl (nperm=1000))
boxplot(disper_ocean)
#sig: p= 0.005

Calcification_Level=betadisper(dist_v, core_metadata2$Calcification_Level, type="centroid")
permutest(Calcification_Level, control=permControl (nperm=1000))
boxplot(Calcification_Level)

disper_cryopreserved=betadisper(dist_v, core_metadata2$cryopreserved, type="centroid")
permutest(disper_cryopreserved, control=permControl (nperm=1000))
boxplot(disper_cryopreserved)

##Checking sifnificance of the model with 3 variables###

rda.all <- rda(core_metadata[,1:457] ~ morphogroup+genus, data=core_metadata, na.action=na.omit)
anova(rda.all, permutations = 1000)
r.adjs <- RsquareAdj(rda.all)$adj.r.squared

summary(rda.all)
vif.cca(rda.all)


rda.numeric <- rda(core_metadata[,1:457] ~ age+Lat+Long+Temperature_maintenance, data=core_metadata, na.action=na.omit)
anova(rda.numeric, permutations = 1000)
r.adjs <- RsquareAdj(rda.numeric)$adj.r.squared

summary(rda.numeric)
vif.cca(rda.numeric)
#should be lower than 10 

##Checking with adonis##

plot(rda.all)
text(rda.all, dis="cn", col="darkblue")
points(rda.all, pch=16, col="black", bg="black", cex=0.8)
text(rda.all,"sites", col="black", cex=0.8, pos=1)
#text("p=0.001", x=80, y=-40)
summary(rda.all)


plot(rda.numeric)
text(rda.numeric, dis="cn", col="darkblue")
points(rda.numeric, pch=16, col="black", bg="black", cex=0.8)
text(rda.numeric,"sites", col="black", cex=0.8, pos=1)
#text("p=0.001", x=80, y=-40)
summary(rda.numeric)
