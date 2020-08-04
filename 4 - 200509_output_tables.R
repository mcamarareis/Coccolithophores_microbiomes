#Script to obtain some information of raw table (number of reads, asvs, etc)

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
library(metacoder); packageVersion("metacoder")

theme_set(theme_bw()) ##define a default theme for ggplot graphics


####Importing tables and performing some summary##
out <- readRDS("200401_out.rds")
seqtab <- readRDS("200406_seqtab.rds")
seqtab_nochimera <- readRDS("200407_seqtab_nochimera.rds")
seqtab_nochimera_newflag <- readRDS("200701_seqtab_nochimera_pooled_newflag.rds")

sum(seqtab)
sum(seqtab_nochimera)
sum(seqtab_nochimera_newflag)

sum(seqtab_nochimera_newflag)/sum(seqtab)

dim(seqtab)
dim(seqtab_nochimera)
dim(seqtab_nochimera_newflag) #more 277 ASVs
dim(seqtab_nochimera_newflag)/dim(seqtab)

table(nchar(getSequences(seqtab)))
table(nchar(getSequences(seqtab_nochimera_newflag)))

###ASVs per sequence legth###

table <- as.data.frame(table(nchar(colnames(seqtab_nochimera_newflag))))
colnames(table) <- c("LENGTH","COUNT")

ggplot(table,aes(x=LENGTH,y=COUNT)) + 
  geom_histogram(stat="identity") + 
  ggtitle("Sequence Lengths by SEQ Count") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=10)) +
  theme(axis.text.y=element_text(size=10))

###Reads per sequence length###

table2 <- tapply(colSums(seqtab_nochimera_newflag), nchar(colnames(seqtab_nochimera_newflag)), sum)
table2 <- data.frame(key=names(table2), value=table2)

colnames(table2) <- c("LENGTH","ABUNDANCE")

ggplot(table2,aes(x=LENGTH,y=ABUNDANCE)) + 
  geom_histogram(stat="identity") + 
  ggtitle("Sequence Lengths by SEQ Abundance") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=10)) +
  theme(axis.text.y=element_text(size=10))

#restricting sequence legth 
tapply(colSums(seqtab_nochimera_newflag), nchar(colnames(seqtab_nochimera_newflag)), sum)
seqlens <- nchar(colnames(seqtab_nochimera_newflag))
seqtab.filt <- seqtab_nochimera_newflag[, (seqlens <= 376) & (seqlens >= 366)]
sum_filt <- tapply(colSums(seqtab.filt), nchar(colnames(seqtab.filt)), sum)

##Obtaining managable names for ASVs

asv_seqs <- colnames(seqtab.filt)
asv_headers <- vector(dim(seqtab.filt)[2], mode="character")

for (i in 1:dim(seqtab.filt)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "200702_ASVs_ehux_restricted_chimera_newflag.fa")

newflag <- read.fasta("200702_ASVs_ehux_restricted_chimera_newflag.fa", forceDNAtolower=FALSE)


# count table:
asv_tab <- t(seqtab.filt)
row.names(asv_tab) <- sub(">", "", asv_headers)
colnames(asv_tab) <- gsub(colnames(asv_tab), pattern = "\\-", replacement = "_") ##replacing "-" by "_" samples names
write.table(asv_tab, "200621_ASVs_counts_ehux_restricted_chimera_newflag.tsv", sep="\t", quote=F, col.names=NA)
asv_tab
sum(asv_tab)
