#### ----- Check Out DEG Results ----- ####

# Load packages 
library(ggplot2)
library(grid)
library(gridBase)
library(VennDiagram)
library(RColorBrewer)
library(tidyverse)
library(ggrepel)

# set working directory 
setwd('~/genetics-tools/thesis/data/')

# Load data 
load('t2dDEG_Drugs.RData')
load('t2dMultiXcan_Drugs.RData')
load('pdDEG_Drugs.RData')
load('pdMultiXcan_Drugs.RData')
load('fda_short.RData')
load('fda_long_ind.RData')
load('fda_long_targ.RData')
load('lincs_spia_out.RData')


# remove binded fda rows to then add fda ind rows 
pdDeg_drugs = pdDeg_drugs[,c(1,2,3,4,5)]
colnames(pdDeg_drugs) = c("Drug","Correlation.Score","Drug.Influenced.Pathway",
                          "Disease.Influenced.Pathway","Affected.Pathway")
pdDeg_drugs = merge(pdDeg_drugs, fda_short, by.x = 'Drug', by.y = 'broad_id')
pdmul_drugs = pdmul_drugs[,c(1,2,3,4,5)]
colnames(pdmul_drugs) = c("Drug","Correlation.Score","Drug.Influenced.Pathway",
                          "Disease.Influenced.Pathway","Affected.Pathway")
pdmul_drugs = merge(pdmul_drugs, fda_short, by.x = 'Drug', by.y = 'broad_id')
t2dDeg_drugs = t2dDeg_drugs[,c(1,2,3,4,5)]
colnames(t2dDeg_drugs) = c("Drug","Correlation.Score","Drug.Influenced.Pathway",
                          "Disease.Influenced.Pathway","Affected.Pathway")
t2dDeg_drugs = merge(t2dDeg_drugs, fda_short, by.x = 'Drug', by.y = 'broad_id')
t2dmul_drugs = t2dmul_drugs[,c(1,2,3,4,5)]
colnames(t2dmul_drugs) = c("Drug","Correlation.Score","Drug.Influenced.Pathway",
                           "Disease.Influenced.Pathway","Affected.Pathway")
t2dmul_drugs = merge(t2dmul_drugs, fda_short, by.x = 'Drug', by.y = 'broad_id')

# subset to get approved drugs 
pdApprDE = subset(pdDeg_drugs, indication == "Parkinson's Disease")
pdApprML = subset(pdmul_drugs, indication == "Parkinson's Disease")
tdApprDE = subset(t2dDeg_drugs, indication == "diabetes mellitus")
tdApprML = subset(t2dmul_drugs, indication == "diabetes mellitus")

#save 
write.table(tdApprML, file = "tdApprML.csv", sep = ",", dec = ".",
            row.names = FALSE, col.names = TRUE)

# Fisher's exact test to test for enrichment 
pdDE_Fish = fisher.test(matrix(c(6,295,17,717),nrow=2,ncol=2))  # pd degs 
pdMUL_Fish = fisher.test(matrix(c(1,300,22,712),nrow=2,ncol=2))  # pd degs 

# plot affected pathways and correlation score for PD DEGS 
model = lm(Correlation.Score~Affected.Pathway, data = pdDeg_drugs)

pdAnnoDE = pdDeg_drugs[,c(1,2,3,4,5)]
colnames(pdAnnoDE) = c("Drug","Correlation.Score","Drug.Influenced.Pathway",
                          "Disease.Influenced.Pathway","Affected.Pathway")
pdAnnoDE = merge(pdAnnoDE, fda_long_ind, by.x = 'Drug', by.y = 'broad_id')
pdAnnoDE = subset(pdAnnoDE, indication == "Parkinson's Disease")
pdAnnoDE = data.frame(
  x = pdAnnoDE$Affected.Pathway,
  y = pdAnnoDE$Correlation.Score, #+ 0.04,
  label = pdAnnoDE$pert_iname
)

pdDeg_drugs %>%
ggplot(aes(x = Affected.Pathway, y = Correlation.Score)) +
  geom_point(alpha = 0.3) + 
  theme(legend.position = "none")

pdDeg_drugs %>%
  ggplot(aes(x = Affected.Pathway, y = Correlation.Score)) +
  geom_point(alpha = 0.4, color = "#3288bd") +
  geom_point(data=pdAnnoDE, 
             aes(x=x,y=y), 
             color="#e41a1c",
             size=2.5) + 
  geom_text_repel(data = pdAnnoDE, aes( x=x, y=y, label=label, fontface = 2),
            color = 'black', 
            size=4) + 
  theme_minimal() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  xlab("% of drug-disease pathway overlap") +
  ylab("Drug-disease correlation") + 
  theme(axis.title=element_text(size=14), 
        axis.text=element_text(size=13)) +
  labs(title = "Parkinson's disease - tissue") +
  theme(plot.title = element_text(face="bold", hjust = 0.5, size = 16)) 

ggsave('pdDE.tiff', width = 6, height = 4, dpi=300)

# plot affected pathways and correlation score for PD MUL
model = lm(Correlation.Score~Affected.Pathway, data = pdmul_drugs)

pdAnnoMUL = pdmul_drugs[,c(1,2,3,4,5)]
colnames(pdAnnoMUL) = c("Drug","Correlation.Score","Drug.Influenced.Pathway",
                       "Disease.Influenced.Pathway","Affected.Pathway")
pdAnnoMUL = merge(pdAnnoMUL, fda_long_ind, by.x = 'Drug', by.y = 'broad_id')
pdAnnoMUL = subset(pdAnnoMUL, indication == "Parkinson's Disease")
pdAnnoMUL = data.frame(
  x = pdAnnoMUL$Affected.Pathway,
  y = pdAnnoMUL$Correlation.Score, 
  label = pdAnnoMUL$pert_iname
)

pdmul_drugs %>%
  ggplot(aes(x = Affected.Pathway, y = Correlation.Score)) +
  geom_point(alpha = 0.3) + 
  theme(legend.position = "none")

pdmul_drugs %>%
  ggplot(aes(x = Affected.Pathway, y = Correlation.Score)) +
  geom_point(alpha = 0.4, color = "#3288bd") +
  geom_point(data=pdAnnoMUL, 
             aes(x=x,y=y), 
             color="#e41a1c",
             size=2.5) + 
  geom_text_repel(data = pdAnnoMUL, aes( x=x, y=y, label=label, fontface = 2),
                  color = 'black', 
                  size=4) + 
  theme_minimal() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  xlab("% of drug-disease pathway overlap") +
  ylab("Drug-disease correlation") + 
  theme(axis.title=element_text(size=14), 
        axis.text=element_text(size=13)) +
  labs(title = "Parkinson's disease - MetaXcan") +
  theme(plot.title = element_text(face="bold", hjust = 0.5, size = 16)) 

ggsave('pdMUL.tiff', width = 6, height = 4, dpi=300)

# plot affected pathways and correlation score for t2d degs
model = lm(Correlation.Score~Affected.Pathway, data = t2dDeg_drugs)

t2dAnnoDE = t2dDeg_drugs[,c(1,2,3,4,5)]
colnames(t2dAnnoDE) = c("Drug","Correlation.Score","Drug.Influenced.Pathway",
                        "Disease.Influenced.Pathway","Affected.Pathway")
t2dAnnoDE = merge(t2dAnnoDE, fda_long_ind, by.x = 'Drug', by.y = 'broad_id')
t2dAnnoDE = subset(t2dAnnoDE, indication == "diabetes mellitus")
t2dAnnoDE = data.frame(
  x = t2dAnnoDE$Affected.Pathway,
  y = t2dAnnoDE$Correlation.Score, 
  label = t2dAnnoDE$pert_iname
)

t2dDeg_drugs %>%
  ggplot(aes(x = Affected.Pathway, y = Correlation.Score)) +
  geom_point(alpha = 0.3) + 
  theme(legend.position = "none")

t2dDeg_drugs %>%
  ggplot(aes(x = Affected.Pathway, y = Correlation.Score)) +
  geom_point(alpha = 0.4, color = "#3288bd") +
  geom_point(data=t2dAnnoDE, 
             aes(x=x,y=y), 
             color="#e41a1c",
             size=2.5) + 
  xlim(0,50) + 
  geom_text_repel(data = t2dAnnoDE, aes( x=x, y=y, label=label, fontface = 2),
                  color = 'black', 
                  size=4) + 
  labs(title = "Type 2 diabetes - tissue") +
  theme(plot.title = element_text(face="bold")) + 
  theme_minimal() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  xlab("% of drug-disease pathway overlap") +
  ylab("Drug-disease correlation") + 
  theme(axis.title=element_text(size=14), 
        axis.text=element_text(size=13)) +
  labs(title = "Type 2 diabetes - tissue") +
  theme(plot.title = element_text(face="bold", hjust = 0.5, size = 16))
ggsave('t2dDE.tiff', width = 6, height = 4, dpi=300)

# plot affected pathways and correlation score for t2d degs
model = lm(Correlation.Score~Affected.Pathway, data = t2dmul_drugs)

t2dAnnoMUL = t2dmul_drugs[,c(1,2,3,4,5)]
colnames(t2dAnnoMUL) = c("Drug","Correlation.Score","Drug.Influenced.Pathway",
                         "Disease.Influenced.Pathway","Affected.Pathway")
t2dAnnoMUL = merge(t2dAnnoMUL, fda_long_ind, by.x = 'Drug', by.y = 'broad_id')
t2dAnnoMUL = subset(t2dAnnoMUL, indication == "diabetes mellitus")
t2dAnnoMUL = data.frame(
  x = t2dAnnoMUL$Affected.Pathway,
  y = t2dAnnoMUL$Correlation.Score, 
  label = t2dAnnoMUL$pert_iname
)

t2dmul_drugs %>%
  ggplot(aes(x = Affected.Pathway, y = Correlation.Score)) +
  geom_point(alpha = 0.3) + 
  theme(legend.position = "none")

t2dmul_drugs %>%
  ggplot(aes(x = Affected.Pathway, y = Correlation.Score)) +
  geom_point(alpha = 0.4, color = "#3288bd") +
  geom_point(data=t2dAnnoMUL, 
             aes(x=x,y=y), 
             color="#e41a1c",
             size=2.5) + 
  xlim(0,50) + 
  geom_text_repel(data = t2dAnnoMUL, aes( x=x, y=y, label=label, fontface = 2),
                  color = 'black', 
                  size=4) + 
  labs(title = "Type 2 diabetes - MetaXcan") +
  theme(plot.title = element_text(face="bold")) + 
  theme_minimal() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  xlab("% of drug-disease pathway overlap") +
  ylab("Drug-disease correlation") + 
  theme(axis.title=element_text(size=14), 
        axis.text=element_text(size=13)) +
  labs(title = "Type 2 diabetes - MetaXcan") +
  theme(plot.title = element_text(face="bold", hjust = 0.5, size = 16))
ggsave('t2dmul.tiff', width = 6, height = 4, dpi=300)


#### ---- END HERE FOR NOW ---- ####
