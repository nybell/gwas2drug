## README for transcriptomic signature matching scripts 
Files and directories listed in the order they need to be run for the analysis. 

/metaxcan-profiles/ - run the scripts in this directory to create the GWAS-based transcriptomic profiles imputed using MetaXcan (2 scripts). 

format-spia-input.R - R script to format all input datato be run in SPIA. Should have 5 files after running this (2 PD transcriptomic signatures, 
                      2 T2D transcriptomic signatures, one file with transcriptomic signatures for each drug). 

spia-lisa.R - R script to run SPIA for the 5 files created above. Needs to be run on a cluster (run using spia-lisa.job)

spia-lisa.job - Job script for running SPIA on lisa

drug-dis-correlations.R - R script for running drug disease correlations for every drug - disease pair. 

tsa-figures.R - R script to create figures. 
