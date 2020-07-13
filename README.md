This repository gathers all the scripts for analysis of E. huxleyi microbiome from raw dereplicated reads, to the final ASV table

# steps_analysis.pdf
Workflow indicating the steps of the analysis performed

# 1 - trimming_primers_cutadapt.qsub

Removing primers using cutadapt. I used anchored primer removal but allowing erros at a rate of 0.1 (it helped to keep more reads - 82% of the total). 

# 2 - 20200114_filtandtrim.R
filterAndTrim step of dada2

# 3 - 20200114_infer_ASVs.R
Learning errors, denoising, chimera removal and taxonomic assignation

# 4 - 200509_output_tables.R
Manipulation of tables and filtering seqs out of legth range

# 5 - 200705_analysis_ehux.R
Performing first analysis of stability, concatenating microbiomes of the same culture and plotting diversity

# 6 - 200712_RDA_ehux.R
Performing PCA, RDA and testing possible drivers of microbiome diversity
