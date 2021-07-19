##### ~~~~~~ Run SPIA for final analysis ~~~~~ #####

# setwd() to data directory 
setwd("~/genetics-tools/thesis/data/")

### ~~~~~~~~~~~~~~~~~~~~~ ###
##### ~~~ LOAD DATA ~~~ #####
### ~~~~~~~~~~~~~~~~~~~~~ ###

# Open Targets data
diaDEGS <- read.csv('targets_DIA.csv', header = T)   # diabetes OT DEGs
parkDEGS <- read.csv('parkinson_degs.csv', header = T)   # parkinson's OT DEGs

# MultiXcan data 
multix_pd <- read.delim("/Users/nyb-macbook/genetics-tools/thesis/data/PD2019-smulti-output.txt", header = TRUE, sep = "\t")    # PD 
multix_t2d <- read.delim("/Users/nyb-macbook/genetics-tools/thesis/data/T2D2018-smulti-output.txt", header = TRUE, sep = "\t")   # T2D


# ----------------------------------------------------- #
#### --- Make MultiXcan outputs SPIA/MAGMA ready --- ####
# ----------------------------------------------------- #

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

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
##### ~~~ FORMAT MULTIXCAN RESULTS FOR SPIA ~~~ #####
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

# clean data using 'make_data_nice.R'

pd_multi_spia_in = list()
for (i in 1 : length(spia_multipd)) {
  pd_multi_spia_in[[i]] = setNames(as.numeric(spia_multipd[[i]][[2]]), as.character(spia_multipd[[i]][[1]]))
  names(pd_multi_spia_in)[[i]] = "Parkinson's Disease"
}

t2d_multi_spia_in = list()
for (i in 1 : length(spia_multit2d)) {
  t2d_multi_spia_in[[i]] = setNames(as.numeric(spia_multit2d[[i]][[2]]), as.character(spia_multit2d[[i]][[1]]))
  names(t2d_multi_spia_in)[[i]] = 'Type II Diabetes Millitus'
}

save(pd_multi_spia_in, file = 'pd_multi_spiainput.RData')
save(t2d_multi_spia_in, file = 't2d_multi_spiainput.RData')

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
##### ~~~ FORMAT PATIENT DEGS FOR SPIA ~~~ #####
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

# need geneID_v97.RData from ...thesis/ps4dr/ps4dr-masters/data/

# assign whatever disease genes to this variable
parkDEGS <-parkDEGS[, -c(1,4,6,7,8,11,12)]
diaDEGS <- diaDEGS[, -c(1,4,6,7,8,11,12)]     # remove extra columns, need ensembl id, gene name,comparison, lfc and pval 

# add efo.id and efo.term for PD
tmp <- parkDEGS
tmp$efo.id <- 'EFO_0002508'                           
tmp$efo.term <- "Parkinson's disease"           
tmp <- tmp[,c("efo.id","efo.term")]
parkDEGS <- cbind(tmp, parkDEGS)
rm(tmp)

# add efo.id and efo.term for T2D
tmp <- diaDEGS
tmp$efo.id <- 'EFO_0001360'                           
tmp$efo.term <- "Type II Diabetes Millitus"           
tmp <- tmp[,c("efo.id","efo.term")]
diaDEGS <- cbind(tmp, diaDEGS)
rm(tmp)

#label col names correctly 
names(parkDEGS) <- c("efo.id","efo.term","ensembl.id","gene.symbol","comparison" ,"lfc","pval")
names(diaDEGS) <- c("efo.id","efo.term","ensembl.id","gene.symbol","comparison" ,"lfc","pval")

# RUN 36-71 TWICE (bad programming but I'm lazy )
DEGs <- data.table(diaDEGS)        # parkDEGS or diaDEGS 

# process positive lfc
lfc_pos = DEGs[lfc >= 0]
#lfc_pos = lfc_pos[order(lfc,decreasing = TRUE), ] #lfc or pvalue we should sort the dataframe?
lfc_pos$pval.pos = abs(log10(lfc_pos$pval)) #created dummy variable for pvalue to order it decreasingly, since there are multiple lfc with same pvalue sometimes.
#lfc_pos = lfc_pos[order(pval), ]
lfc_pos = lfc_pos[order(lfc, pval.pos, decreasing = TRUE),]
lfc_pos = lfc_pos[! duplicated(lfc_pos[, c('efo.id', 'ensembl.id')]),]
#lfc_pos.efo = split(lfc_pos, lfc_pos$efo.id)
lfc_pos$pval.pos = NULL

# process negative lfc
lfc_neg = DEGs[lfc < 0]
#lfc_neg = lfc_neg[order(lfc), ] #lfc or pvalue we should sort the dataframe?
lfc_neg = lfc_neg[order(pval, lfc),]
lfc_neg = lfc_neg[! duplicated(lfc_neg[, c('efo.id', 'ensembl.id')]),]
#lfc_neg.efo = split(lfc_neg, lfc_neg$efo.id)

# combine both lfc
lfc_comb = rbind(lfc_pos, lfc_neg) #combine positive and negative lfc change back to a single data frame
lfc_comb = lfc_comb[order(pval),]
lfc_comb = lfc_comb[! duplicated(lfc_comb[, c('efo.id', 'ensembl.id')]),]
lfc_comb$pval = NULL
rm(lfc_neg, lfc_pos, DEGs)

# add ensembl and entrez ids 
lfc_comb = merge(lfc_comb, gene_id, by = "ensembl.id")
lfc_comb = lfc_comb[! duplicated(lfc_comb[, c('efo.id', 'ENTREZ')]),]

# drop the comparison and gene symbol columns 
lfc_comb <- lfc_comb[,-c(4,5)]

# split into list 
lfc_t2d = split(lfc_comb, lfc_comb$efo.term) # change to lfc_t2d for diabetes 

# should extract efo.term, log fold change, and entrez gene id's and format for SPIA 
pd_deg_spia_in = list()
for (i in 1 : length(lfc_pd)) {
  pd_deg_spia_in[[i]] = setNames(as.numeric(lfc_pd[[i]][[4]]), as.character(lfc_pd[[i]][[5]]))
  names(pd_deg_spia_in)[[i]] = names(lfc_pd)[[i]]
}

t2d_deg_spia_in = list()
for (i in 1 : length(lfc_t2d)) {
  t2d_deg_spia_in[[i]] = setNames(as.numeric(lfc_t2d[[i]][[4]]), as.character(lfc_t2d[[i]][[5]]))
  names(t2d_deg_spia_in)[[i]] = names(lfc_t2d)[[i]]
}

save(pd_deg_spia_in, file = 'pd_deg_spiainput.RData')
save(t2d_deg_spia_in, file = 't2d_deg_spiainput.RData')


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
##### ~~~ Format LINCS Signatures For SPIA ~~~ #####
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

# keep only top and bottom 15 expressed genes

topgene_sigs = list()
for (i in seq(exemp_mat)) {
  tmp = exemp_mat[[i]]
  lincs_id <- names(exemp_mat)[i]
  tmp_top = slice_max(tmp, zmat, n = 15)
  tmp_bot = slice_min(tmp, zmat, n = 15)
  tmp_all = rbind(tmp_top, tmp_bot)
  topgene_sigs[i] = list(tmp_all)
  names(topgene_sigs)[i] <- lincs_id
}

rm(tmp, tmp_all, tmp_bot, tmp_top)

exmp_sigs = topgene_sigs

#loop through every drug in fda short and row bind z-score

lincs_spia_input = list()
for (i in fda_short$broad_id) {
  lincs_spia_input[[i]] = setNames(as.numeric(exmp_sigs[[i]][['zmat']]), as.character(exmp_sigs[[i]][['entrez']]))
}

names(lincs_spia_input) = names(exmp_sigs)

# preprocess by removing drugs genes with no overlap with entrez all 
#entrez_uni = data.frame(entrez_all)
#lincs_spia_inputV2 = list()
#for (i in seq(lincs_spia_input)) {
#  tmp1 = data.frame(names(lincs_spia_input[[i]]))
#  tmp2 = data.frame(lincs_spia_input[[i]])
#  tmp3 = cbind(tmp1,tmp2)
#  colnames(tmp3) = c('entrez', 'zmat')
#  foo = merge(entrez_uni, tmp3, by.x = 'entrez_all', by.y = 'entrez')
#  lincs_spia_inputV2[i] = list(foo)
#  names(lincs_spia_inputV2)[i] = names(lincs_spia_input)[i]
#}

#lincs_spia_inputV3 = list()
#for (i in seq(lincs_spia_inputV2)) {
#  lincs_spia_inputV3[[i]] = setNames(as.numeric(lincs_spia_inputV2[[i]][['zmat']]), as.character(lincs_spia_inputV2[[i]][['entrez_all']]))
#}


# add pert names to lincs_spia_input
#names(lincs_spia_input) = fda_short$pert_iname
#names(lincs_spia_inputV3) = names(lincs_spia_inputV2)

save(lincs_spia_input, file = 'lincs_spia_input.RData')
save(lincs_spia_inputV3, file = 'lincs_spia_inputV3.RData')

#### --- FINISHED: FINAL NOTES --- ####

# Now take the five SPIA input files and plug them into the "spia-lisa.R" and use the 
# spia-lisa.job script to run SPIA on LISA
