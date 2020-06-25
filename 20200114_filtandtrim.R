library(dada2); packageVersion("dada2")

TRUNCLEN_1=215
TRUNCLEN_2=190

# Input fastq primer trimmed directory

path <-"/projet/umr7144/mapp/mcamaradosreis/ehuxleyi_microbiome/archives/all_data/primer_trimmed"
##########

# List fastqs 
fastqFs <- sort(list.files(path, pattern="_R1_trimmed.fastq.gz",  full.names = TRUE))
fastqRs <- sort(list.files(path, pattern="_R2_trimmed.fastq.gz",  full.names = TRUE))

# Filtered fastqs name
filtFs <- gsub("primer_trimmed", "quality_filtered",fastqFs)
filtRs <- gsub("primer_trimmed", "quality_filtered",fastqRs)

# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
out <- filterAndTrim(fastqFs, filtFs, fastqRs, filtRs, truncLen=c(TRUNCLEN_1,TRUNCLEN_2),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)

saveRDS(out, "200407_out.rds")
