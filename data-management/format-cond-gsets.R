#### ---- Format gene-sets for conditional analysis ---- ####

# load libraries
library(dplyr)
library(biomaRt)
library('org.Hs.eg.db')

# set working directory 
setwd('~/genetics-tools/thesis/data/')

# Load data 
load('fda_short.RData')
load('fda_long_ind.RData')
load('fda_long_targ.RData')


#### ----  make drug gene sets from drug targets (Gaspar method) ---- ####

symbols <- as.vector(fda_long_targ$target)                        # create vector of gene symbols
entrez_ids <- mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')    # create equivalent vector with entrez ids 
ensembl_ids <- mapIds(org.Hs.eg.db, symbols, 'ENSEMBL', 'SYMBOL')    # create equivalent with ensembl ids

# add ensembl and entrez ids to fda_long_targets
fda_long_targ$ensembl <- ensembl_ids
fda_long_targ$entrez <- entrez_ids

# format to save as gene sets 
gene_set_entrez <- fda_long_targ[, c('pert_iname', 'entrez')]

# collapse rows by lincs_id & combine data (1 row = 1 drug) for patient DEGS
entrez_short <- gene_set_entrez%>%group_by(pert_iname)%>%summarise(entrez = paste(entrez, collapse=","))%>%mutate(entrez = strsplit(as.vector(entrez), ","))

# druggable genes 
dGenes = data.frame(unique(fda_long_targ$entrez))
colnames(dGenes) = 'X'
dGenes = dGenes%>%summarise(X = paste(X, collapse=","))%>%mutate(X = strsplit(as.vector(X), ","))
cond_drug = cbind('cond_drug', dGenes)  
colnames(cond_drug) = c('pert_iname', 'entrez')
cond_gene_sets = rbind(entrez_short, cond_drug)

# make pretty 
cond_gene_sets$entrez = gsub("^c\\(", "", cond_gene_sets$entrez)
cond_gene_sets$entrez = gsub("\\)", "", cond_gene_sets$entrez)
cond_gene_sets$entrez = gsub("\"", "", cond_gene_sets$entrez)
cond_gene_sets$entrez = trimws(cond_gene_sets$entrez, which = "both")
write.table(cond_gene_sets, file = 'cond_gene_sets.txt', sep = "\t", col.names = TRUE, row.names = FALSE)

# 
dGenes = data.frame(unique(fda_long_targ$entrez))
dGenes = rbind('cond_drug', dGenes)

moaSets = rbind.fill(moa_uniq, data.frame(t(dGenes$unique.fda_long_targ.entrez.)))
indSets = rbind.fill(ind_uniq, data.frame(t(dGenes$unique.fda_long_targ.entrez.)))

write.table(moaSets, file = 'moa_cond_gsets.txt', sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(indSets, file = 'ind_cond_gsets.txt', sep = "\t", col.names = FALSE, row.names = FALSE)

### --- END --- ###

