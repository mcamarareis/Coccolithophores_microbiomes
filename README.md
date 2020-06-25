This repository gathers all the scripts for analysis of E. huxleyi microbiome from raw dereplicated reads, to the final ASV table

# steps_analysis.pdf
Workflow indicating the steps of the analysis performed

# trimming_primers_cutadapt.qsub

Removing primers using cutadapt. I used anchored primer removal but allowing erros at a rate of 0.1 (it helped to keep more reads - 82% of the total). 

# 20200114_filtandtrim.R
filterAndTrim step of dada2

# 20200114_infer_ASVs.R
Learning errors, denoising, chimera removal and taxonomic assignation

# 200509_pre_processing_ehux.R
Manipulation of tables and some analysis in R
