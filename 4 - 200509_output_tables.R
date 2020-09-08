#' ---
#' title: "Obtaining statistics of the rawdata and saving tables for downstream analyzes"
#' author: "Mariana Câmara dos Reis"
#' date: "September 8th, 2020"
#' output:
#'    html_document:
#'      toc: true
#'      highlight: haddock
#' ---

#' This script is the first for manipulating the ASVs table.
#' First we obtain the number of reads and ASVs before and after chimera removal.
#' We check as well the legths distribution of the ASVs and we filter those out of the range of legth of the primers used.
#' Next we export the **.fasta** and **.csv** files that will be used for downstream analyzes.


#+ echo = FALSE 
setwd("/Users/marianacamaradosreis/Documents/PhD/Chapter1/Final_results/analysis_R/200331_microbiome_analysis")

#' # Load the libraries
#+ messages = FALSE
library(ggplot2); packageVersion("ggplot2")
library(dplyr); packageVersion("dplyr")
library(tidyverse); packageVersion("tidyverse")
library(dada2) ;packageVersion("dada2")
library(readr);packageVersion("readr")
library(ape);packageVersion("ape")
theme_set(theme_bw()) ##define a default theme for ggplot graphics

#' # Importing the output files from dada2

seqtab <- readRDS("200406_seqtab.rds") #the ASV table
seqtab_nochimera_newflag <- readRDS("200701_seqtab_nochimera_pooled_newflag.rds") #the ASV table after removing chimera

#' # Checking number of reads 
#' Before chimera removal

sum(seqtab)

#' After chimera removal

sum(seqtab_nochimera_newflag)

#' Checking % of reads kept

sum(seqtab_nochimera_newflag)/sum(seqtab)
#' 90% of the reads were kept

#' # Checking number of ASVs
#' Before chimera removal

dim(seqtab)

#' After chimera removal

dim(seqtab_nochimera_newflag)

#' Checking % of ASVs remaining

dim(seqtab_nochimera_newflag)/dim(seqtab)
#' The majority of the ASVs was removed as chimera

#' # Checking sequences length distribution
#' Before chimera removal

table(nchar(getSequences(seqtab)))

#' After chimera removal

table(nchar(getSequences(seqtab_nochimera_newflag)))

#' # Plotting the distribution of sequence lengths
 
#' ## Number of ASVs

table <- as.data.frame(table(nchar(colnames(seqtab_nochimera_newflag))))
colnames(table) <- c("LENGTH","COUNT")
#+ fig.width=9, fig.height=6
ggplot(table,aes(x=LENGTH,y=COUNT)) + 
  geom_histogram(stat="identity") + 
  ggtitle("Sequence Lengths by SEQ Count") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=10)) +
  theme(axis.text.y=element_text(size=10))

#' ## Number of reads

table2 <- tapply(colSums(seqtab_nochimera_newflag), nchar(colnames(seqtab_nochimera_newflag)), sum)
table2 <- data.frame(key=names(table2), value=table2)

colnames(table2) <- c("LENGTH","ABUNDANCE")
#+ fig.width=9, fig.height=6
ggplot(table2,aes(x=LENGTH,y=ABUNDANCE)) + 
  geom_histogram(stat="identity") + 
  ggtitle("Sequence Lengths by SEQ Abundance") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=10)) +
  theme(axis.text.y=element_text(size=10))

#' Since most of the reads and ASVs are restricted between 366 and 376 bp I'll restric sequences to this legth.
#' In general sequences out of the range of sequence length of the primers (in our case the primers give 373 bp) are related to contamination, unespecific primming or eukaryotic sequence (here chloroplasts and mithocondria).
#' Since they correspond to only a very small proportion of the total of reads (less than 1%), I'll remove it.

#' # Filtering out the sequences larger than 376 and smaller than 366 bp

tapply(colSums(seqtab_nochimera_newflag), nchar(colnames(seqtab_nochimera_newflag)), sum)
seqlens <- nchar(colnames(seqtab_nochimera_newflag))
seqtab.filt <- seqtab_nochimera_newflag[, (seqlens <= 376) & (seqlens >= 366)]
sum_filt <- tapply(colSums(seqtab.filt), nchar(colnames(seqtab.filt)), sum)

#' # Obtaining managable names for ASVs

asv_seqs <- colnames(seqtab.filt)
asv_headers <- vector(dim(seqtab.filt)[2], mode="character")

for (i in 1:dim(seqtab.filt)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

#' # Exporting the fasta file with the sequences kept

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "200702_ASVs_ehux_restricted_chimera_newflag.fa")

#' # Exporting the count table

asv_tab <- t(seqtab.filt)
row.names(asv_tab) <- sub(">", "", asv_headers)
colnames(asv_tab) <- gsub(colnames(asv_tab), pattern = "\\-", replacement = "_") ##replacing "-" by "_" samples names

write.table(asv_tab, "200621_ASVs_counts_ehux_restricted_chimera_newflag.tsv", sep="\t", quote=F, col.names=NA)

#' # Checking the table

head(asv_tab)

#' # Checking the number of reads kept after filtering the sequences out of the range of legth

sum(asv_tab)/sum(seqtab_nochimera_newflag)

#' 99% of the reads were kept

sessionInfo()
