This repository gathers all the scripts for analysis of coccolithophores microbiomes from raw demultiplexed reads, to the final ASV table and statistical analysis.

# steps_analysis.pdf
Workflow indicating the steps of the analysis performed

# 1 - trimming_primers_cutadapt.qsub

Removing primers using cutadapt. 

# 2 - 
Dada2 workflow
# 3 - 20200114_infer_ASVs.R
Learning errors, denoising, chimera removal and taxonomic assignation

# 4 - output_tables.R
Manipulation of tables and filtering seqs out of legth range
html_report: 

# 5 - 
Performing first analysis of stability, concatenating microbiomes of the same culture and plotting diversity

# 6 - 
Performing PCA, RDA and testing possible drivers of microbiome diversity
