##Script to compare remove chloroplasts, mitochondria, and to study stability of microbiome over time##

##Heading from dada2 to phyloseq##
setwd("/Users/marianacamaradosreis/Documents/PhD/Chapter1/Final_results/analysis_R/200331_microbiome_analysis")

##to install new necessary packages:  BiocManager::install("name of package")"


library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(data.table); packageVersion("data.table")
library(vegan); packageVersion("vegan")
library(pheatmap); packageVersion("pheatmap")
library(dplyr); packageVersion("dplyr")
library(tidyverse); packageVersion("tidyverse")
library(RColorBrewer); packageVersion("RColorBrewer")
library(dada2) ;packageVersion("dada2")
library(readr);packageVersion("readr")
library(ape);packageVersion("ape")
detach(package:plyr) ####phyr changes the behavior of grouping in dyplr. So If I use ggbiplot, I have to open a new session after.   
library(dplyr)
library(metacoder); packageVersion("metacoder")

theme_set(theme_bw()) ##define a default theme for ggplot graphics

###Manipulating tha data using tidyverse####

dt <- read_tsv("200621_ASVs_counts_ehux_restricted_chimera_newflag.tsv")
print(dt)
dim(dt)
sum(dt[2:223])

tax_data <- read_tsv("200705_ASVs_TAX_ehux_restricted_chimera_newflag_IDTaxa_50.tsv") %>%
  replace(., is.na(.), "unknown")
print(tax_data)
dim(tax_data)

metadata <- read_tsv("/Users/marianacamaradosreis/Documents/PhD/Chapter3/Final_results/200505_metadata_simplified.txt")  
print(metadata)
dim(metadata)

colnames(dt)[2:223] == metadata$id
dt$X1 == tax_data$ASVId

dt2 <- left_join(dt, tax_data,
                      by = c("X1" = "ASVId"))

tail(colnames(dt2), n = 20)
print(dt2$genus)


##1) Filtering out Chloroplasts, Mitochondria and Eukaryotes###
####empty cells need to be replaced with "unknown" otherwise they are discarded

dt3 <- dt2 %>% 
  filter(., domain == 'Bacteria' & class != 'Cyanobacteriia' & family != 'Mitochondria')
dim(dt3)
sum(dt3[2:223])

#1515 5259947

##2) Obtaining some stats (number of reads per sample, etc) before and after removing chloroplasts

before <- dt %>% 
  gather(var, value, -X1) %>% 
  mutate_at(vars(value), list(as.numeric)) %>%
  mutate(Total = sum(value)) %>%
  group_by(var) %>% 
  summarise(Reads_per_sample = sum(value))


after <-  dt3[,1:223] %>% 
  gather(var, value, -X1) %>% 
  mutate_at(vars(value), list(as.numeric)) %>%
  mutate(Total = sum(value)) %>%
  group_by(var) %>% 
  summarise(Reads_per_sample_nochloroplast = sum(value))

tracking_reads <- inner_join(before, after, by = c("var" = "var"))
dim(tracking_reads)

write.table(tracking_reads, "200703_reads_per_sample.txt", sep="\t")
##3) Checking ordination with all samples

dtp <- dt3[,1:223] %>%
  rownames_to_column %>% 
  gather(var, value, -rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value)

colnames(dtp) <- dtp[1,]
dtp <- dtp[-1,,drop=F]

dtpm <- inner_join(dtp, metadata, by = c("X1" = "id")) %>% 
  filter(keep_ehux_with_out == "yes")
tail(colnames(dtpm), n=10)


###For this I remove ASVs with less than 0.001% of the total of reads. Which in this case is: 52.6

metadata <- read_tsv("/Users/marianacamaradosreis/Documents/PhD/Chapter3/Final_results/200505_metadata_simplified.txt")  %>%
  filter(keep_ehux_with_out== "yes")
print(metadata)

dtpm2 <- dtpm[,1:1515] %>% 
  gather(var, value, -X1) %>% 
  mutate_at(vars(value), list(as.numeric)) %>% 
  group_by(var) %>% 
  mutate(ASV_abund = sum(value)) %>% 
  filter(ASV_abund > 52)

min(dtpm2$ASV_abund)

dtpm_hell <- dtpm2 %>%
  select(X1, value, var) %>% 
  pivot_wider(names_from = var, values_from = value) %>% 
  replace(., is.na(.), 0) %>% column_to_rownames(., var = "X1") %>% 
  decostand(., method="hellinger") 

dim(dtpm_hell)
#483
dist_b <- vegdist(dtpm_hell, "bray")
hc <- hclust(dist_b, "complete")
hc$labels <- metadata$rcc



##4) Doing the same, but after removing outliers

dtp <- dt3[,1:223] %>%
  rownames_to_column %>% 
  gather(var, value, -rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value)

colnames(dtp) <- dtp[1,]
dtp <- dtp[-1,,drop=F]

dtpm <- inner_join(dtp, metadata, by = c("X1" = "id")) %>% 
  filter(keep_ehux == "yes")
tail(colnames(dtpm), n=10)

dim(dtpm)


###For this I remove ASVs with less than 0.001% of the total of reads. Which in this case is: 52.6

metadata <- read_tsv("/Users/marianacamaradosreis/Documents/PhD/Chapter3/Final_results/200505_metadata_simplified.txt")  %>%
  filter(keep_ehux == "yes")
print(metadata)

dtpm2 <- dtpm[,1:1515] %>% 
  gather(var, value, -X1) %>% 
  mutate_at(vars(value), list(as.numeric)) %>% 
  mutate(Total = sum(value)) %>%
  group_by(var) %>% 
  mutate(ASV_abund = sum(value)) %>% 
  filter(ASV_abund > 52)

min(dtpm2$ASV_abund)
dim(dtpm2)

dtpm_hell <- dtpm2 %>%
  select(X1, value, var) %>% 
  pivot_wider(names_from = var, values_from = value) %>% 
  replace(., is.na(.), 0) %>% column_to_rownames(., var = "X1") %>% 
  decostand(., method="hellinger") 

dim(dtpm_hell)

#466

dist_b <- vegdist(dtpm_hell, "bray")
hc <- hclust(dist_b, "complete")
hc$labels <- metadata$rcc

plot(hc, cex=0.7)

##5) Exporting number of reads after abundance filtering###
dtpm_abund_filt <- dtpm2 %>%
  select(X1, value, var) %>% 
  pivot_wider(names_from = var, values_from = value) %>% 
  replace(., is.na(.), 0) %>% column_to_rownames(., var = "X1") 

abund_filt <- dtpm_abund_filt %>% 
  rownames_to_column(., "X1") %>%
  gather(var, value, -X1) %>% 
  mutate_at(vars(value), list(as.numeric)) %>%
  mutate(Total = sum(value)) %>%
  group_by(X1) %>% 
  summarise(Reads_per_sample_abund_filt = sum(value))

write.table(abund_filt, "200703_reads_per_sample_abund_filt.txt", sep="\t")

##6) Plotting the mean distance for each culture

metadata <- read_tsv("/Users/marianacamaradosreis/Documents/PhD/Chapter3/Final_results/200505_metadata_simplified.txt")  %>%
  filter(keep_ehux == "yes")
print(metadata)

colnames(as.matrix(dist_b))==metadata$id

meanDist <- meandist(dist_b, metadata$rcc_2)
meanDist
summary(meanDist)


##Taking the diagonal of the matrix ##

diag=diag(x = meanDist, names=TRUE)
diag <- as.matrix(diag)

df <- data.frame(diag) %>% rownames_to_column(., "strain")

df

ggplot(df) +
  labs(y = "Mean Bray-Curtis Dissimilarity", cex=16) +
  geom_bar( aes(x=strain, y=diag), stat="identity", fill="gray47", alpha=0.7) +
 #geom_errorbar( aes(x=strain, ymin=diag-diag_sd, ymax=diag+diag_sd), width=0.4, colour="black", alpha=0.9) +
  theme_classic(base_size = 16) + theme(axis.text.x = element_text(face="bold", size=12, angle=90), axis.text.y = element_text(face="bold", size=12))

##7) Merging samples based on average of present ASVs##
#important to use na.rm=FALSE, it will not do mean if of the replicates has one na
###Using th abundance filtered table###
dim(dtpm)

dtpm2 <- dtpm[,1:1515] %>% 
  gather(var, value, -X1) %>% 
  mutate_at(vars(value), list(as.numeric)) %>% 
  mutate(Total = sum(value)) %>%
  group_by(var) %>% 
  mutate(ASV_abund = sum(value)) %>% 
  filter(ASV_abund > 52)

dtpm3 <- dtpm2 %>%
  select(X1, value, var) %>% 
  pivot_wider(names_from = var, values_from = value) %>% 
  replace(., is.na(.), 0)


dtpm_forgrouping <- inner_join(dtpm3, metadata, by = c("X1" = "id")) %>% 
  filter(keep_ehux == "yes") %>% na_if(., "0") %>% mutate_at(vars(starts_with("ASV_")),list(as.numeric))


group <- group_by(dtpm_forgrouping, rcc)  %>% summarise(across(starts_with("ASV_"), .fns = mean, na.rm=FALSE)) %>%
  replace(., is.na(.), 0) %>% 
  mutate_at(vars(starts_with("ASV_")),list(as.numeric)) %>%
  column_to_rownames(., "rcc") %>%
  .[, which(colSums(.) > 0)]

write.table(group, "200705_core_ehux_merged_newflag.txt", sep="\t")

##7) Importing grouped table back for next analyses

group <- read_tsv("200705_core_ehux_merged_newflag.txt")
print(group)
dim(group)
sum(group[,2:458])
#1846514

###Preparing data to plot with taxonomy##
##What I need, data filtered by abundance but without zero#

toplot <- group %>% gather(var, value, -X1) %>% 
  group_by(var) %>%
  mutate(Total = sum(value)) 

min(toplot$Total)

###Diversity heatree
library(ggplot2)
library(viridis)
library(hrbrthemes)  
library(ggthemes)

toplot <- toplot %>% group_by(X1) %>% 
  mutate(Total_sample = sum(value)) %>% 
  mutate(Relative_abund=(value/Total_sample)*100)

tax_data <- read_tsv("200705_ASVs_TAX_ehux_restricted_chimera_newflag_IDTaxa_50.tsv") %>%
  replace(., is.na(.), "unclassified")
print(tax_data)

metadata_conc <- read_tsv("200712_metadata_strains.txt")
print(metadata_conc)

toplot_wide <- toplot %>% select(X1, value, var) %>% 
  pivot_wider(names_from = X1, values_from = value) %>% 
  replace(., is.na(.), 0) #%>% column_to_rownames(., var = "var")
dim(toplot_wide)

toplot_tax <- left_join(toplot_wide, tax_data,
                        by = c("var" = "ASVId")) %>%
  column_to_rownames(., var = "var")
dim(toplot_tax)
sample_data$sample == colnames(toplot_tax)[1:76]


tail(colnames(toplot_tax), n=20)
toplot_tax[,77:82]
obj <- parse_tax_data(toplot_tax,
                      class_cols = 77:82,
                      named_by_rank = TRUE) 

obj$data$tax_abund <- calc_taxon_abund(obj, data = "tax_data", 
                                       cols = metadata$sample)


print(obj$data$tax_abund)


##Plotting taxonomic data -->heattree

set.seed(1) # This makes the plot appear the same each time it is run 
taxa_patterns_to_remove <- paste0("^", "unclassified", "$")
obj %>% 
  mutate_obs("tax_abund", abundance = (rowSums(obj$data$tax_abund[,2:77])/1846324)) %>%
  filter_taxa(taxon_ranks == "family", supertaxa = TRUE) %>% # subset to the order rank
  filter_taxa(! Reduce(`|`, lapply(taxa_patterns_to_remove, grepl, x = taxon_names))) %>%
  #filter_taxa(n_obs > 1) %>%
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color = abundance,
            node_color_axis_label = "Relative abundance",
            node_size_axis_label = "ASV count",
            initial_layout = "re", layout = "da")


##8) Obtaining some stats of the main groups###

abund_tax_family <- inner_join(toplot, tax_data, by=c("var" = "ASVId")) %>%
  mutate_at(vars(value), list(as.numeric)) %>%
  group_by(family) %>% 
  summarise(family_abund = sum(value)) %>%
  mutate(Total = sum(family_abund)) %>%
  mutate(family_relative_abund = (family_abund/Total)*100) %>% 
  arrange(., desc(as.numeric(family_relative_abund)))

write.table(abund_tax_family, "200705_summary_family.tsv", sep="\t", quote=F, col.names=NA)

abund_tax_phylum <- inner_join(toplot, tax_data, by=c("var" = "ASVId")) %>%
  mutate_at(vars(value), list(as.numeric)) %>%
  group_by(phylum) %>% 
  summarise(phylum_abund = sum(value)) %>%
  mutate(Total = sum(phylum_abund)) %>%
  mutate(phylum_relative_abund = (phylum_abund/Total)*100) %>% 
  arrange(., desc(as.numeric(phylum_relative_abund)))

write.table(abund_tax_phylum, "200705_summary_phylum.tsv", sep="\t", quote=F, col.names=NA)

abund_mainASVS <- inner_join(toplot, tax_data, by=c("var" = "ASVId")) %>%
  mutate_at(vars(value), list(as.numeric)) %>%
  group_by(var) %>% 
  summarise(asv_abund = sum(value)) %>%
  mutate(Total = sum(asv_abund)) %>%
  mutate(asv_relative_abund = (asv_abund/Total)*100) %>% 
  arrange(., desc(as.numeric(asv_relative_abund)))

write.table(abund_mainASVS, "200705_summary_asvs.tsv", sep="\t", quote=F, col.names=NA)


