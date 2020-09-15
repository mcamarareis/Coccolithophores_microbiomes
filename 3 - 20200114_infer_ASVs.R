library(dada2); packageVersion("dada2")

# File parsing 

filtpath <- "/projet/umr7144/mapp/mcamaradosreis/ehuxleyi_microbiome/archives/all_data/quality_filtered"

filtFs <- list.files(filtpath, pattern="R1_trimmed.fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpath, pattern="R2_trimmed.fastq.gz", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) 
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) 

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names

set.seed(100)

##Learning errors##
###################
errF <- learnErrors(filtFs, nbases = 2e+08, randomize=TRUE)
errR <- learnErrors(filtRs, nbases = 2e+08, randomize=TRUE)

##ASVs inference##
##################

dadaFs <- dada(filtFs, err=errF, pool=TRUE)
dadaRs <- dada(filtRs, err=errR, pool=TRUE)

####Merging pairs###
####################

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

##exporting all intermediary files
saveRDS(errF, "200407_errF.rds")
saveRDS(errR, "200407_errR.rds")
saveRDS(mergers,"200407_mergers.rds")
saveRDS(dadaFs, "200407_ddF.rds")
saveRDS(dadaRs, "200407_ddR.rds")

# Construct sequence table 
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "200406_seqtab.rds")

sum(seqtab)

#Chimera removal#

seqtab_nochimera <- removeBimeraDenovo(seqtab, method="pooled", minFoldParentOverAbundance=8, verbose=TRUE)
sum(seqtab_nochimera)
saveRDS(seqtab_nochimera, "200701_seqtab_nochimera_pooled_newflag.rds")


#Reads trhough the pipeline
out <- readRDS("200407_out.rds")

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab_nochimera))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Assign taxonomy
taxa <- assignTaxonomy(seqtab_nochimera, "/projet/umr7144/mapp/mcamaradosreis/ehuxleyi_microbiome/archives/all_data/quality_filtered/silva_nr_v132_train_set.fa")








