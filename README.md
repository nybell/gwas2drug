## gwas2drug ##

This GitHub repository contains all the data, results, and code used during my master's thesis. For this project we compared 2 drug repurposing methods, transcriptomic signature analysis and drug gene-set analysis. 

# Overview of directories

data-management: Contains R code to format a variety of the data used to match the input requirements for various software tools used (e.g., MetaXcan, MAGMA)

data: contains both raw input data and result data files for both drug repurposing analysis. Several files are missing due to size. 

drug-gene-set-analysis: contains the scripts used to run the drug gene-set analysis in MAGMA. 

lincs-code: code to extract the LINCS L1000 transcriptomic signatures using their cmapR package. This needs to be run on a cluster (LISA). 

transcriptomic-signature-analysis: contains all the code to run the transcriptomic signature matching analysis. 

NOTE: many sub-directory README files are unfinished, and will be updated shortly. 
