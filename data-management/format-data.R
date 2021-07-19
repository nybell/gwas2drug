##### ----- Format data for different purposes ------ #####



#### set working directory ####
setwd("/Users/nyb-macbook/genetics-tools/thesis/data/")

#### load packages ####
library(stringr)
library(data.table)
library(dplyr)
library(tidyr)
library(pracma)
library(ggplot2)
library(grid)
library(gridBase)
library(biomaRt)
library('org.Hs.eg.db')

# load data
meta_clue <- read.delim("drug-meta-clue.txt", header = TRUE, sep = "\t")       # load clue repo meta 
sample_clue <- read.delim("drug-sample-clue.txt", header = TRUE, sep = "\t")   # load clue sample data
load("/Users/nyb-macbook/genetics-tools/thesis/data/fda_sigs2020.RData")       # load drug DEGs
multix_pd <- read.delim("/Users/nyb-macbook/genetics-tools/thesis/data/PD2019-smulti-output.txt", header = TRUE, sep = "\t")
multix_t2d <- read.delim("/Users/nyb-macbook/genetics-tools/thesis/data/T2D2018-smulti-output.txt", header = TRUE, sep = "\t")

##### ~~~ Create data set with all drug data in one place ('fda_drugs') ~~~ #####
# merge pert names from clue_sample with pert names from clue_drug
fda_drugs <- merge(sample_clue[, c("broad_id", "pert_iname", "pubchem_cid"),],meta_clue, by = "pert_iname")   
fda_drugs$broad_id = (data.frame(gsub('.{9}$', '', fda_drugs$broad_id)))
fda_drugs = data.table(fda_drugs)
colnames(fda_drugs) <- c("pert_iname", "broad_id","pubchem_cid","clinical_phase","moa","target","disease_area","indication")  
fda_drugs = fda_drugs[!duplicated(fda_drugs), ]

# subset fda_drugs to only contain drugs in our sample
our_lincs <- data.frame(names(exemp_mat))
colnames(our_lincs) <- 'broad_id'
fda <- merge(our_lincs, fda_drugs, by = 'broad_id')

# create fda short (one drug per row)
fda_tmp = subset(fda, indication != "")
fda_short = subset(fda_tmp, target != "")
rm(fda_tmp, our_lincs)

# remove duplicates from fda_short
ind = !duplicated(fda_short$pert_iname)
fda_short = fda_short[ind,]

# unnest clue indication to make fda_long_ind (each row is a unique drug-indication pair)
fda_long_ind = fda_short %>%
  mutate(indication=strsplit(indication, "\\|")) %>% 
  unnest(indication)

# unnest gene targets to make fda_long_targ (each row = unique drug-gene target pair)
fda_long_targ = fda_short %>%
  mutate(target=strsplit(target, "\\|")) %>% 
  unnest(target)

#### -----  make drug gene sets from drug targets (drug gene-set method) ------ ###

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

# make pretty 
entrez_short$entrez = gsub("^c\\(", "", entrez_short$entrez)
entrez_short$entrez = gsub("\\)", "", entrez_short$entrez)
entrez_short$entrez = gsub("\"", "", entrez_short$entrez)
entrez_short$entrez = trimws(entrez_short$entrez, which = "both")

# merge with fda short (so entrez ids are added)
fda_short = merge(fda_short, entrez_short, by = 'pert_iname')

# save all of the above
write.table(entrez_short, file = 'entrez_genesets.txt', sep = "\t", col.names = FALSE, row.names = FALSE)    # save entrez gene-sets for Gaspar method
write.table(fda_short, file = 'fda_short.txt', sep = "\t", col.names = TRUE, row.names = FALSE)    # save fda_short.txt
save(fda_short, file = 'fda_short.RData')             # save fda_short.RData
write.table(fda_long_ind, file = 'fda_long_ind.txt', sep = "\t", col.names = TRUE, row.names = FALSE)   # save fda_long_ind.txt
save(fda_long_ind, file = 'fda_long_ind.RData')       # save fda_long_ind.RData
write.table(fda_long_targ, file = 'fda_long_targ.txt', sep = "\t", col.names = TRUE, row.names = FALSE) # save fda_long_targ.txt
save(fda_long_targ, file = 'fda_long_targ.RData')     # save fda_long_targ.RData

#### ----- make multixcan outputs SPIA/MAGMA ready ----- ####

# remove .* from each ensembl gene name for pd multi 
multix_pd$gene <- gsub("\\..*","",multix_pd$gene)
multix_pd = data.table(multix_pd)

# remove .* from each ensembl gene name for pd multi 
multix_t2d$gene <- gsub("\\..*","",multix_t2d$gene)
multix_t2d = data.table(multix_t2d)

# map ensembl ids to entrez ids
pdmulti_entrez_ids <- mapIds(org.Hs.eg.db, multix_pd$gene, 'ENTREZID', 'ENSEMBL')    # create equivalent vector with entrez ids 
t2dmulti_entrez_ids <- mapIds(org.Hs.eg.db, multix_t2d$gene, 'ENTREZID', 'ENSEMBL')    # create equivalent vector with entrez ids 

# cbind to output 
multix_pd <- cbind(pdmulti_entrez_ids, multix_pd)
multix_t2d <- cbind(t2dmulti_entrez_ids, multix_t2d)

# check number of na values for mapped genes
print(sum(is.na(multix_pd$pdmulti_entrez_ids)))
print(sum(is.na(multix_t2d$t2dmulti_entrez_ids)))

# remove incomplete rows (NA entrez ids)_pd = complete.cases(multix_pd)
multix_pd = multix_pd[complete.cases(multix_pd)]
multix_t2d = multix_t2d[complete.cases(multix_t2d)]

#bonferroni correct
multix_pd$bonf = p.adjust(multix_pd$pvalue, method = 'bonferroni', n = length(multix_pd$pvalue))
multix_t2d$bonf = p.adjust(multix_t2d$pvalue, method = 'bonferroni', n = length(multix_t2d$pvalue))

# remove nonsignificant rows
multiPD_bonf = subset(multix_pd, bonf < 0.05)
multiT2D_bonf = subset(multix_t2d, bonf < 0.05)

# remove duplicates from fda_short
ind = !duplicated(multiT2D_bonf$t2dmulti_entrez_ids)
multiT2D_bonf = multiT2D_bonf[ind,]

# make spia inputs 
spia_multipd <- as.vector(multiPD_bonf[, c('pdmulti_entrez_ids', 'z_mean')])
spia_multit2d <- as.vector(multiT2D_bonf[, c('t2dmulti_entrez_ids', 'z_mean')])
spia_multipd <- list(spia_multipd)
spia_multit2d <- list(spia_multit2d)

#### ----- format magma lincs gene sets (unused) ----- #####

# making exemp mats into data frame
entrez_lincs <- exemp_mat[["BRD-K76022557"]][["entrez"]]
gene_sym_lincs <- exemp_mat[["BRD-K76022557"]][["gene_symbol"]]
ens_lincs <- exemp_mat[["BRD-K76022557"]][["ensembl_id"]]

df_lincs <- rbind(entrez_lincs, gene_sym_lincs, ens_lincs)

#loop through every drug in fda short and row bind z-score 
for (i in fda_short$broad_id) {
  tmp = exemp_mat[[i]][['zmat']]
  df_lincs <- rbind(df_lincs, tmp)
}
rm(tmp,i)

# transpose to 1 column = transcriptomic sig for drug, then add drug names 
tmp = data.frame(t(df_lincs[4:1038,]))
colnames(tmp) = fda_short$pert_iname

# bind to gene labels 
lincs_entrez_gsets <- cbind(entrez_lincs, tmp)
# lincs_sym_gsets <- cbind(gene_sym_lincs, tmp)
# lincs_ens_gsets <- cbind(ens_lincs, tmp)
rm(tmp, df_lincs)

# write lincs drug gene sets to text file 
write.table(lincs_entrez_gsets, file = 'lincs_entrez_gsets.txt', sep = "\t", col.names = TRUE, row.names = FALSE)

#### ----- Format gene sets based on grouping drug by mechanism of action ----- ####

# unnest gene targets to make fda_long_targ (each row = unique drug/moa pair)
fda_long_moa = fda_long_targ %>%
  mutate(moa=strsplit(moa, "\\|")) %>% 
  unnest(moa)

# make data frame of only moa and entrez
entrez_moa = fda_long_moa[,c('moa', 'entrez')]

# split by moa 
moa_list = split(entrez_moa, entrez_moa$moa)

# loop through and remove duplicates in gene rows
for (i in seq(moa_list)) {
  if (i == 1) {
    drug = names(moa_list)[i]
    genes = t(unique(moa_list[[i]][["entrez"]]))
    moa_uniq = cbind(drug, genes)
    length(moa_uniq) = 100
  } else {
    drug = names(moa_list)[i]
    genes = t(unique(moa_list[[i]][["entrez"]]))
    tmp = cbind(drug, genes)
    length(tmp) = 100
    moa_uniq = rbind(moa_uniq, tmp)
  }
}
#make into df
moa_uniq = data.frame(moa_uniq)

# add - 
moa_uniq$X1 = gsub(" ","-", moa_uniq$X1)

# save
write.table(moa_uniq, file = 'moa_uniq_gsets.txt', sep = "\t", col.names = FALSE, row.names = FALSE)

#### ----- Format gene sets based on grouping drug by indication ----- ####

# unnest gene targets to make fda_long_targ (each row = unique drug/moa pair)
fda_tar_ind = fda_long_targ %>%
  mutate(indication=strsplit(indication, "\\|")) %>% 
  unnest(indication)

# make data frame of only moa and entrez
entrez_ind = fda_tar_ind[,c('indication', 'entrez')]

# split by moa 
ind_list = split(entrez_ind, entrez_ind$indication)

# loop through and remove duplicates in gene rows
for (i in seq(ind_list)) {
  if (i == 1) {
    drug = names(ind_list)[i]
    genes = t(unique(ind_list[[i]][["entrez"]]))
    ind_uniq = cbind(drug, genes)
    length(ind_uniq) = 250
  } else {
    drug = names(ind_list)[i]
    genes = t(unique(ind_list[[i]][["entrez"]]))
    tmp = cbind(drug, genes)
    length(tmp) = 250
    ind_uniq = rbind(ind_uniq, tmp)
  }
}
#make into df
ind_uniq = data.frame(ind_uniq)

# add - 
ind_uniq$X1 = gsub(" ","-", ind_uniq$X1)

# save
write.table(ind_uniq, file = 'ind_uniq_gsets.txt', sep = "\t", col.names = FALSE, row.names = FALSE)


#### ---- END ---- ####
