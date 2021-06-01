#' ---
#' title: "Checking final number of reads and outputing tables"
#' author: "Mariana Câmara dos Reis"
#' date: "`r Sys.Date()`"
#' output:
#'   html_document:
#'     fig_retina: NULL
#'     highlight: haddock
#'     toc: true
#'     toc_float:
#'       collapsed: false
#'       smooth_scroll: false
#'     code_folding: hide
#'     code_download: true
#' ---
#' 
#+ echo=F
knitr::opts_chunk$set(warning=FALSE,message=F)

#' This script is the first for manipulating the ASVs table.
#' First we obtain the number of reads and ASVs before and after chimera removal.
#' We check as well the legths distribution of the ASVs and we filter those out of the range of legth of the primers used.
#' Next we export the **.fasta** and **.csv** files that will be used for downstream analyzes.


#Load the libraries
#+ messages = FALSE
library(ggplot2); packageVersion("ggplot2")
library(dplyr); packageVersion("dplyr")
library(tidyverse); packageVersion("tidyverse")
library(dada2) ;packageVersion("dada2")
library(readr);packageVersion("readr")
library(ape);packageVersion("ape")
library(data.table);packageVersion("data.table")
library(knitr);packageVersion("knitr")
library(DECIPHER); packageVersion("DECIPHER")
theme_set(theme_bw()) ##define a default theme for ggplot graphics

#' ## Plotting error models and checking convergence


errF=readRDS("my_favorite_directory/201206_errF.rds")
errR=readRDS("my_favorite_directory/201206_errR.rds")

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

dada2:::checkConvergence(errF)
dada2:::checkConvergence(errR)

#' Convergence after 5 rounds for forward and after 4 for reverse. Error models are good.

seqtab <- readRDS("my_favorite_directory/201206_seqtab.rds") #the ASV table
seqtab_nochimera <- readRDS("my_favorite_directory/201206_seqtab_nochimera.rds") #the ASV table after removing chimera

#' ## Checking number of reads{.tabset}
#' Before chimera removal
sum(seqtab)
#' After chimera removal
sum(seqtab_nochimera)
#' Checking % of reads kept
sum(seqtab_nochimera)/sum(seqtab)
#' 90% of the reads were kept

#' ## Checking number of ASVs
#' Before chimera removal
dim(seqtab)
#' After chimera removal
dim(seqtab_nochimera)
#' Checking % of ASVs remaining
dim(seqtab_nochimera)/dim(seqtab)
#' The majority (70%) of the ASVs was removed as chimera

#' ## Checking sequences length distribution
#' Before chimera removal
table(nchar(getSequences(seqtab)))

#' After chimera removal
table(nchar(getSequences(seqtab_nochimera)))

#' ## Plotting the distribution of sequence lengths  {.tabset}
 
#' ### Number of ASVs

table <- as.data.frame(table(nchar(colnames(seqtab_nochimera))))
colnames(table) <- c("LENGTH","COUNT")
#+ fig.width=9, fig.height=6
ggplot(table,aes(x=LENGTH,y=COUNT)) + 
  geom_histogram(stat="identity") + 
  ggtitle("Sequence Lengths by SEQ Count") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=10)) +
  theme(axis.text.y=element_text(size=10))

#' ### Number of reads

table2 <- tapply(colSums(seqtab_nochimera), nchar(colnames(seqtab_nochimera)), sum)
table2 <- data.frame(key=names(table2), value=table2)

colnames(table2) <- c("LENGTH","ABUNDANCE")
#+ fig.width=9, fig.height=6
ggplot(table2,aes(x=LENGTH,y=ABUNDANCE)) + 
  geom_histogram(stat="identity") + 
  ggtitle("Sequence Lengths by SEQ Abundance") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=10)) +
  theme(axis.text.y=element_text(size=10))

#' ## Selecting ASVs out of reads length to check their taxonomy

tapply(colSums(seqtab_nochimera), nchar(colnames(seqtab_nochimera)), sum)
seqlens <- nchar(colnames(seqtab_nochimera))
larger <- seqtab_nochimera[, (seqlens > 376)]
smaller <- seqtab_nochimera[, (seqlens < 366)]

all(rownames(larger)==rownames(smaller))

removed <- cbind(larger, smaller)
dim(removed)

#' 136 ASVs are out of the range of length
removed_seqs <- colnames(removed)
removed_headers <- vector(dim(removed)[2], mode="character")

for (i in 1:dim(removed)[2]) {
  removed_headers[i] <- paste(">ASV", i, sep="_")
}

removed_fasta <- c(rbind(removed_headers, removed_seqs))
write(removed_fasta, "tmp/ASVs_removed.fa")

#' # Exporting the count table

removed_tab <- t(removed)
row.names(removed_tab) <- sub(">", "", removed_headers)
colnames(removed_tab) <- gsub(colnames(removed_tab), pattern = "\\-", replacement = "_") ##replacing "-" by "_" samples names
asvs_removed_abund <- data.frame(rowSums(removed_tab)) %>%
  data.table(., keep.rownames = T)
write.table(removed_tab, "tmp/ASVs_counts_removed.tsv", sep="\t", quote=F, col.names=NA)
kable(asvs_removed_abund)

#' ## Assgining taxonomy of these ASVs with vsearch
#' 
 
fas <- "tmp/ASVs_removed.fa"
reference_515F_926R <- "my_favorite_directory/silva_v138/SILVA_138_SSURef_NR99_tax_silva_515F_926R.fasta"
tax_removed <- "tmp/tax_ehux_removed.tsv"

system(paste0("/Users/marianacamaradosreis/Downloads/vsearch-2.15.0-macos-x86_64/bin/vsearch ",
              "--usearch_global ",  fas, 
              " --db ",              reference_515F_926R,
              " --id 0.80 ",
              " --blast6out ",          tax_removed))

tax_removed <- fread(paste0(tax_removed))
tax_removed_abund <- left_join(asvs_removed_abund, tax_removed, by=c("rn"="V1"))
write.csv(tax_removed_abund, "tmp/tax_ehux_removed.tsv")

#' Usually the sequences out of the range of length of the primers appear by unespecific priming, contamination of eukaryotic sequences. 
#' By checking the otu table of these sequences I can see that: The most abundant sequences are not assigned to any organism, which seems to be bad sequences (even using blast).
#' One sequence that is assigned to bacteria appear randomly in few cultures (but not in all the replicates of the same culture), so it will be removed by the filter of prevalence. 
#' The others are very low abundance ASVs (in addition of not being classified) which would be lost by abundance threshold, so I 'll remove them now.

#' ## Filtering out the sequences larger than 376 and smaller than 366 bp

tapply(colSums(seqtab_nochimera), nchar(colnames(seqtab_nochimera)), sum)
seqlens <- nchar(colnames(seqtab_nochimera))
seqtab.filt <- seqtab_nochimera[, (seqlens <= 376) & (seqlens >= 366)]
sum_filt <- tapply(colSums(seqtab.filt), nchar(colnames(seqtab.filt)), sum)

#' ## Obtaining managable names for ASVs

asv_seqs <- colnames(seqtab.filt)
asv_headers <- vector(dim(seqtab.filt)[2], mode="character")

for (i in 1:dim(seqtab.filt)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}


#' ## Exporting the fasta file with the sequences kept

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "my_favorite_directory/ASVs_coccolithophores_restricted_chimera.fa")

#' ## Exporting the count table

asv_tab <- t(seqtab.filt)
row.names(asv_tab) <- sub(">", "", asv_headers)
colnames(asv_tab) <- sapply(strsplit(colnames(asv_tab), "_"), `[`, 1)
colnames(asv_tab) <- gsub(colnames(asv_tab), pattern = "\\-", replacement = "_") ##replacing "-" by "_" samples names
write.table(asv_tab, "my_favorite_directory/ASVs_counts_coccolithophores_restricted_chimera.tsv", sep="\t", quote=F, col.names=NA)


#' ## Filtering tax table and exporting with ASVs names
tax_data <- fread("my_favorite_directory/201206_ASVs_tax_coccolithophores_microbiomes.tsv") %>%
  replace(., is.na(.), "unclassified")


tax_filtered <- tax_data[tax_data$V1 %in% colnames(seqtab.filt),]
dim(tax_filtered)
asv_headers <- vector(dim(tax_filtered)[1], mode="character")

for (i in 1:dim(tax_filtered)[1]) {
  asv_headers[i] <- paste("ASV", i, sep="_")
}

tax_filtered$ASVId <- asv_headers
write.table(tax_filtered, "my_favorite_directory/ASVs_TAX_coccolithophores_restricted.tsv", sep="\t", quote=F, col.names=NA)


#' ## Checking the table

kable(head(asv_tab))

#' ## Checking the number of reads kept after filtering the sequences out of the range of legth

sum(asv_tab)/sum(seqtab_nochimera)

#' 99% of the reads were kept

sessionInfo()
