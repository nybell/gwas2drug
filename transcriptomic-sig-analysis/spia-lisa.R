#### ---- Run SPIA for each disease profile ---- ####

##### Load packages #####
library(cowplot)
library(data.table)
library(doParallel)
library(dplyr)
library(gridExtra)
library(Hmisc)
library(httr)
library(jsonlite)
library(pROC)
library(purrr)
library(stringr)
library(tools)
library(tidyr)
library(BiocParallel)
library(biomaRt)
library(graphite)
library(SPIA)
library(org.Hs.eg.db)

# -------------: SPIA Functions :------------- #

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
}


spia_fun = function(x){                            # not working, need to change variable names
  spia_kegg_drug <- foreach(i = seq(x)) %dopar% {
    foreach(j = seq(x[[i]])) %dopar% {
      drug = x[[i]][[j]]
      tmp=spia(de = drug, all = entrez_all, data.dir = file.path(dataFolder,"spia_input/real_kegg/"), organism = "hsa")
    }
  }
}

name_fun = function(spia_kegg_drug){              # not working, need to change names variable names
  for (disease in 1:length(spia_kegg_drug)) {
    names(spia_kegg_drug)[[disease]] = names(lfc_test)[[disease]]
    for (drug in 1:length(spia_kegg_drug[[disease]])) {
      names(spia_kegg_drug[[disease]])[[drug]] = names(lfc_test[[disease]])[[drug]]
    }
  }
  spia_kegg_drug
}

# ------------- define dataFolder pathway ------------- #

#dataFolder <- '/home/nbell/ps4dr/ps4dr_master/data/'
dataFolder <- '~/genetics-tools/thesis/ps4dr/ps4dr_master/data/'

# ------------- load 'entrez_all' ------------- #

load("geneID_v97.RData")
entrez_all = unique(gene_id$ENTREZ)
print("entrez all loaded")

# ------------- setwd() & load data ------------- # 

setwd('/home/nbell/ps4dr/spia/input/')

load('pd_deg_spiainput.RData')         # PD patient DEGs 
load('t2d_deg_spiainput.RData')        # type 2 patient DEGs 
load('pd_multi_spiainput.RData')       # PD multixcan 
load('t2d_multi_spiainput.RData')      # type 2 multixcan 
load('lincs_spia_input.RData')         # drug data 

# ------------- create ref gene sets ------------- #

lincsAll = as.integer(names(lincs_spia_input[[1]]))
pdMultiALL = as.integer(multix_pd$pdmulti_entrez_ids)
t2dMultiALL = as.integer(multix_t2d$t2dmulti_entrez_ids)

# ------------- DISEASE KEGG SPIA ------------- #

# NOTE: if getting error message from SPIA make sure that there are no duplicated genes in the input

# run on PD DEGs
pdDEGS_spia_out = list()
for (i in 1 : length(pd_deg_spia_in)) {
  pdDEGS_spia_out[[i]] = spia(de = pd_deg_spia_in[[i]], all = entrez_all, data.dir = file.path(dataFolder,"spia_input/real_kegg/"), organism = "hsa")
}

names(pdDEGS_spia_out) = names(pd_deg_spia_in)
save(pdDEGS_spia_out, file = "spia_out_pddegs.RData")
print('SPIA PD DEGs Done')

# run on t2d DEGs
t2dDEGS_spia_out = list()
for (i in 1 : length(t2d_deg_spia_in)) {
  t2dDEGS_spia_out[[i]] = spia(de = t2d_deg_spia_in[[i]], all = entrez_all, data.dir = file.path(dataFolder,"spia_input/real_kegg/"), organism = "hsa")
}

names(t2dDEGS_spia_out) = names(t2d_deg_spia_in)
save(t2dDEGS_spia_out, file = "spia_out_t2ddegs.RData")
print('SPIA T2D DEGs Done')

# run on pd multixcan
pdmul_spia_out = list()
for (i in 1 : length(pd_multi_spia_in)) {
  pdmul_spia_out[[i]] = spia(de = pd_multi_spia_in[[i]], all = pdMultiALL, data.dir = file.path(dataFolder,"spia_input/real_kegg/"), organism = "hsa")
}

names(pdmul_spia_out) = names(pd_multi_spia_in)
save(pdmul_spia_out, file = "spia_out_pdmul.RData")
print('SPIA PD MultiXcan Done')

# run on t2d multixcan
t2dmul_spia_out = list()
for (i in 1 : length(t2d_multi_spia_in)) {
  t2dmul_spia_out[[i]] = spia(de = t2d_multi_spia_in[[i]], all = t2dMultiALL, data.dir = file.path(dataFolder,"spia_input/real_kegg/"), organism = "hsa")
}

names(t2dmul_spia_out) = names(t2d_multi_spia_in)
save(t2dmul_spia_out, file = "spia_out_t2dmul.RData")
print('SPIA T2D MultiXcan Done')

print("Disease SPIA finished")

# -------------: Drug SPIA - KEGG :------------- #

lincs_spia_out = list()
for (i in seq(lincs_spia_input)) {
  lincs_spia_out[[i]] = spia(de = lincs_spia_input[[i]], all = lincsAll, data.dir = file.path(dataFolder,"spia_input/real_kegg/"), organism = "hsa")
  print(i)
}
names(lincs_spia_out) = names(lincs_spia_input)

save(lincs_spia_out, file = "lincs_spia_out.RData")

print("Drug SPIA finished")


#### ---- END ---- ####
