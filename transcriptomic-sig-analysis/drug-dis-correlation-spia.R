#### ----- SPIA Drug Disease Correlations ----- ####

# NOTE: need to run script separately for each phenotype. Change lines 39-40 and 67-70 accordingly. 

#### ----- Setup ---- ####

# load packages
library(data.table)
library(tidyr)
library(dplyr)
library(RecordLinkage)
library(ggplot2)
library(gridExtra)
library(Hmisc)
library(purrr)
library(tools)
library(plotly)
library(cowplot)
library(qqplotr)
library(pROC)
library(rlist)

# set working directory 
setwd('~/genetics-tools/thesis/data/')

# load data 
load('lincs_spia_out.RData')
load('spia_out_pdmul.RData')
load('spia_out_t2dmul.RData')
load('spia_out_t2ddegs.RData')
load('spia_out_pddegs.RData')
load('fda_short.RData')
load('fda_long_ind.RData')

#### ----- Clean Data ----- ####

#' load Kegg Drug SPIA file
spia_drug = list(lincs_spiaOut)
names(spia_drug) = "Parkinson's disease"   # might need to be "Parkinson's disease" depending on input file 
#names(spia_drug) = "Type II Diabetes Millitus"

spia_drug = Filter(function(x) ! is.null(x), spia_drug) #delete empty df from list

#~~~~Remove any drug pathway with p.value (pNDE) >= 0.05 ~~~#

for (i in seq_along(spia_drug)) {
  spia_drug[[i]] = lapply(spia_drug[[i]], function(x) x[x$pNDE <= 0.05,])
  spia_drug[[i]] = spia_drug[[i]][lapply(spia_drug[[i]], length) > 1]
  spia_drug[[i]] = Filter(function(x) ! dim(x)[1] == 0, spia_drug[[i]])
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drug_path_temp = vector('list', length(spia_drug)) # create list of lists
names(drug_path_temp) = names(spia_drug)
for (i in seq_along(spia_drug)) {
  for (j in seq_along(spia_drug[[i]])) {
    drug_path_temp[[i]][[j]] = spia_drug[[i]][[j]][, c(1, 11)] # use 2 for ID
    names(drug_path_temp[[i]])[[j]] = names(spia_drug[[i]])[[j]]
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~: load KEGG Disease SPIA results :~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

spia_disease = pdDEGS_spia_out
#spia_disease = pdmul_spia_out
#spia_disease = t2dDEGS_spia_out
#spia_disease = t2dmul_spia_out

# Create disease paths with only pathway name and their activity status
# Remove any other diseases which are not in drug_path_temp
# so, both drug_path_temp & dis_path are equivalent

dis_path = vector('list', length(drug_path_temp)) # create list of lists
names(dis_path) = names(drug_path_temp)

for (i in seq_along(drug_path_temp)) {
  for (j in seq_along(spia_disease)) {
    if (names(drug_path_temp)[[i]] == names(spia_disease)[j]) {
      dis_path[[i]] = spia_disease[[j]][, c(1, 11)]
    }
  }
}

#' Now, remove those diseases from drug_path_temp which are not in dis_path

drug_path = vector('list', length(dis_path)) # create list of lists
names(drug_path) = names(dis_path)

for (i in 1 : length(drug_path_temp)) {
  if (names(drug_path_temp)[[i]] %in% names(spia_disease)) {
    drug_path[[i]] = drug_path_temp[[i]]
    names(drug_path)[[i]] = names(drug_path_temp)[[i]]
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~: Calculate Correlation-Score :~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drug_dis_path = vector('list', length(drug_path)) # create list of lists
names(drug_dis_path) = names(drug_path)
drug_correlation = vector('list', length(drug_path)) # create list of lists
names(drug_correlation) = names(drug_path)

# assign +1 and -1 for values "activated" and "inhibited"
for (i in seq_along(drug_path)) {
  for (j in seq_along(drug_path[[i]])) {
    drug_dis_path[[i]][[j]] = merge(dis_path[[i]], drug_path[[i]][[j]], by = "Name")
    names(drug_dis_path[[i]])[[j]] = names(drug_path[[i]])[[j]]
    names(drug_dis_path[[i]][[j]]) = c("Pathways", "Disease.Influence", "Drug.Influence")
    drug_dis_path[[i]][[j]]$Disease.Influence = ifelse(drug_dis_path[[i]][[j]]$Disease.Influence == "Activated", 1, - 1)
    drug_dis_path[[i]][[j]]$Drug.Influence = ifelse(drug_dis_path[[i]][[j]]$Drug.Influence == "Activated", 1, - 1)
  }
}

# run the correlations to get p-value (excludes ~1/2 of drugs in data)
#leftout_drugs = list()
#for (i in seq_along(drug_path)) {
#  for (j in seq_along(drug_path[[i]])) {
#    if (length(drug_dis_path[[i]][[j]]$Drug.Influence) < 3) {
#      leftout_drugs[[j]] = drug_dis_path[[i]][[j]]
#      names(leftout_drugs)[[j]] = names(drug_dis_path[[i]])[[j]]
#      next
#    } else {
#      drug_correlation[[i]][[j]] = cor.test(drug_dis_path[[i]][[j]]$Drug.Influence, drug_dis_path[[i]][[j]]$Disease.Influence)[["estimate"]][["cor"]]   
#      names(drug_correlation[[i]])[[j]] = names(drug_dis_path[[i]])[[j]]
#      drug_correlation[[i]][[j]] = as.data.frame(drug_correlation[[i]][[j]])
#      names(drug_correlation[[i]][[j]]) = "Correlation.Score"
#      drug_correlation[[i]][[j]]$"Correlation.pvalue" = cor.test(drug_dis_path[[i]][[j]]$Drug.Influence, drug_dis_path[[i]][[j]]$Disease.Influence)$p.value
#      drug_correlation[[i]][[j]]$"DrugPathway" = length(drug_dis_path[[i]][[j]]$Drug.Influence)
#      drug_correlation[[i]][[j]]$"DiseasePathway" = length(dis_path[[i]]$Name)
#      drug_correlation[[i]][[j]]$"affectedPathway" = round((drug_correlation[[i]][[j]]$"DrugPathway" / drug_correlation[[i]][[j]]$"DiseasePathway") * 100, 2)
#      drug_correlation[[i]][[j]]$"Disease" = names(dis_path)[[i]]
#    }
#  }
#}

# run using PS4DR method (don't use p-value)
for (i in seq_along(drug_path)) {
  for (j in seq_along(drug_path[[i]])) {
    drug_correlation[[i]][[j]] = cor(drug_dis_path[[i]][[j]]$Drug.Influence, drug_dis_path[[i]][[j]]$Disease.Influence)   
    names(drug_correlation[[i]])[[j]] = names(drug_dis_path[[i]])[[j]]
    drug_correlation[[i]][[j]] = as.data.frame(drug_correlation[[i]][[j]])
    names(drug_correlation[[i]][[j]]) = "Correlation.Score"
    #drug_correlation[[i]][[j]]$"Pathways" = paste(drug_dis_path[[i]][[j]][["Pathways"]], sep = ",")     # want pathway data 
    drug_correlation[[i]][[j]]$"DrugPathway" = length(drug_dis_path[[i]][[j]]$Drug.Influence)
    drug_correlation[[i]][[j]]$"DiseasePathway" = length(dis_path[[i]]$Name)
    drug_correlation[[i]][[j]]$"affectedPathway" = round((drug_correlation[[i]][[j]]$"DrugPathway" / drug_correlation[[i]][[j]]$"DiseasePathway") * 100, 2)
    drug_correlation[[i]][[j]]$"Disease" = names(dis_path)[[i]]
  }
}

#' create a single data.frame from all drugs in a disease

for (i in seq_along(drug_correlation)) {
  drug_correlation[[i]] = do.call(rbind, drug_correlation[[i]])
  drug_correlation[[i]]$Drug = rownames(drug_correlation[[i]])
  drug_correlation[[i]] = drug_correlation[[i]][, c(6, 1, 2, 3, 4, 5)]
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

drug_correlation = do.call(rbind, drug_correlation) # list to df
drug_correlation = filter(drug_correlation, ! is.na(Correlation.Score)) # remove all rows with correlation score NA
drug_correlation = split(drug_correlation, drug_correlation$Disease) # df to list
drug_correlation = lapply(drug_correlation, data.table) # make all df to data table

drug_correlation = Filter(function(x) dim(x)[1] >= 1, drug_correlation) # remove empty disease lists

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#' filter out drugs from each disease with correlationscore greater than -0.1
drug_shortlist = drug_correlation
drug_shortlist = lapply(drug_shortlist, data.table)

#' remove drugs with less than the given thresholds
for (i in seq_along(drug_shortlist)) {
  drug_shortlist[[i]] = drug_shortlist[[i]][drug_shortlist[[i]]$Correlation.Score < 0]
}

#' Order all lists in drug_shortlist
for (i in seq_along(drug_shortlist)) {
  drug_shortlist[[i]] = drug_shortlist[[i]][order(drug_shortlist[[i]][['Correlation.Score']],-drug_shortlist[[i]][['affectedPathway']]),]
}

drug_shortlist_df = do.call(rbind, drug_shortlist)
drug_shortlist_df = drug_shortlist_df[,c(1,2,3,4,5)]
names(drug_shortlist_df) = c("Drug","Correlation.Score","Drug.Influenced.Pathway","Disease.Influenced.Pathway","Affected.Pathway(%)")

# correct for multiple tests
#drug_shortlist_df$bonf = p.adjust(drug_shortlist_df$P.value, method = 'bonferroni', n = length(drug_shortlist_df$P.value))
#drug_shortlist_df = subset(drug_shortlist_df, P.value < 0.05)
#drugDF_bonf = subset(drug_shortlist_df, bonf < 0.05)


### merge with fda short and save ###

#NOTE: Need to run script separately for each phenotype. Un-comment code below corresponding to the data used.

# PD MultiXcan Results 
#pdmul_drugs = merge(drug_shortlist_df, fda_short, by.x = 'Drug', by.y = 'broad_id')
#save(pdmul_drugs, file = 'pdMultiXcan_Drugs.RData')

# PD DEG Results 
#pdDeg_drugs = merge(drug_shortlist_df, fda_short, by.x = 'Drug', by.y = 'broad_id')
#pdDeg_drugs = merge(drugDF_bonf, fda_short, by.x = 'Drug', by.y = 'broad_id')     # bonf corrected 
#save(pdDeg_drugs, file = 'pdDEG_Drugs.RData')

# Type 2 Diabetes MultiXcan Results 
#t2dmul_drugs = merge(drug_shortlist_df, fda_short, by.x = 'Drug', by.y = 'broad_id')
#save(t2dmul_drugs, file = 't2dMultiXcan_Drugs.RData')

# Type 2 Diabetes DEGs Results 
#t2dDeg_drugs = merge(drug_shortlist_df, fda_short, by.x = 'Drug', by.y = 'broad_id')
#save(t2dDeg_drugs, file = 't2dDEG_Drugs.RData')





### --- END --- ###
