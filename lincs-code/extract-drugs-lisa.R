##### ~~~ LINCs L1000 Pre-Processing ~~~ #####

# set working directory 
setwd("/home/nbell/lincs-data/")

# load packages
library(pracma)
library(BiocManager)
library(cmapR)
library(dplyr)
library(tidyr)
library(data.table)

# create pathway variables, make sure setwd() is correct 
ds_path <- "/home/nbell/lincs-data/level5_beta_trt_cp_n720216x12328.gctx"
siginfo_path <- "siginfo_beta.txt"
geneinfo_path <- "geneinfo_beta.txt"
cellinfo_path <- "cellinfo_beta.txt"
compoundinfo_path <- "compoundinfo_beta.txt"

# read signature annotations (corresponding to columns of level 5 matrix)
siginfo <- data.table::fread(siginfo_path)
cellinfo <- data.table::fread(cellinfo_path)
compounndinfo <- data.table::fread(compoundinfo_path)
geneinfo <- data.table::fread(geneinfo_path)

# load in clue drug repurposing meta data 
clue_drug <- read.delim("drug-meta-clue.txt", header = TRUE, sep = "\t")
clue_sample <- read.delim("drug-sample-clue.txt", header = TRUE, sep = "\t")

# reformat broad id's in clue_sample to be 'core ids'
brd_ids <- unique(data.frame(gsub('.{9}$', '', clue_sample$broad_id)))
brd_ids = data.table(brd_ids) 
colnames(brd_ids) <- 'core_ids'

# merge pert names from clue_sample with pert names from clue_drug
fda_drugs <- merge(clue_sample[, c("broad_id", "pert_iname", "pubchem_cid"),],clue_drug, by = "pert_iname")   
fda_drugs$broad_id = (data.frame(gsub('.{9}$', '', fda_drugs$broad_id)))
fda_drugs = data.table(fda_drugs)
colnames(fda_drugs) <- c("pert_iname", "broad_id","pubchem_cid","clinical_phase","moa","target","disease_area","indication")  
fda_drugs = fda_drugs[!duplicated(fda_drugs), ]


##### ~~~~~ get drug signatures ~~~~~ #####

# loop through every drug we have have clinical information on & extract signature info 
pre_drug_sigs = list()
for (i in seq(fda_drugs$broad_id)) {
  tmp = siginfo[pert_id == fda_drugs$broad_id[i]]
  pre_drug_sigs[i] = list(tmp)
  names(pre_drug_sigs)[i] <- fda_drugs$broad_id[i]
  rm(tmp)
}

# remove emtpy lists
drug_sigs = list()
for (i in seq(pre_drug_sigs)) { 
  if (isempty(pre_drug_sigs[[i]][["sig_id"]])) {
    next 
  } else {
    if (sum(pre_drug_sigs[[i]][["is_ncs_exemplar"]]) >= 1) {
      tmp <- subset(pre_drug_sigs[[i]], is_ncs_exemplar == 1)
      drug_sigs[i] <- list(tmp)
      names(drug_sigs)[i] <- names(pre_drug_sigs[i])
    } else if (pre_drug_sigs[[i]][["qc_pass"]] >= 1) {
      tmp <- subset(pre_drug_sigs[[i]], qc_pass == 1)
      drug_sigs[i] <- list(tmp)
      names(drug_sigs)[i] <- names(pre_drug_sigs[i])
    } else {
      drug_sigs[i] <- pre_drug_sigs[i]
      names(drug_sigs)[i] <- names(pre_drug_sigs[i])
    }
  }
}

# remove null lists 
drug_sigs = drug_sigs[lengths(drug_sigs) != 0]

# loop through and extract data for each 
ex_list = list()
for (i in seq(drug_sigs)) { 
  tmp = drug_sigs[[i]]
  ex_list[i] <- parse_gctx(ds_path,
                           cid = tmp$sig_id)
  rm(tmp)
}

# add names 
names(ex_list) = names(drug_sigs)

# merge all gene IDs and labels into one data.frame
rid = data.frame(ex_list[[1]]@rid)
colnames(rid) = 'entrez'
gene_data <- merge(rid, geneinfo[, c("gene_id","gene_symbol", "ensembl_id")], by.x = 'entrez', by.y = 'gene_id')

# loop through ex_list, extract matrix & rowMeans to get average gene expression across signatures
exemp_mat = list()
for (i in seq(ex_list)){
  zmat = mat(ex_list[[i]])
  zmat = data.frame(zmat)
  zmat = rowMeans(zmat)
  zmat = cbind(gene_data, zmat)
  exemp_mat[i] = list(data.frame(zmat))
  rm(zmat)
}

# add names
names(exemp_mat) = names(drug_sigs)

#save 
save(exemp_mat, file = "fda_sigs2020.RData")


