setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("4_Chemostat.R")
library(EnhancedVolcano)

level <-"Species"
#############################################Section one : litter 123###################################################################################
ps <- subset_samples(PSB, Sample_time_point %in% c("Chemostat4"))

ps <- tax_glom(ps, level , NArm = FALSE) #select a level to compare
#draw heatmap 

#Include deseq2 p.adj  taxonomy for rows
#load two list of deseq2
tab_1_all <- read.table("Bacteriome/stat_result/Chemostat/Deseq/Lactose_HMO_chemostat4.tsv", sep = "\t",header = T, row.names = 1)

###Volcano plot#####################################
#Load packages needed for making volcano plot
#Cut off padj=0.05, log2FoldChange=0.6
p_1<-  EnhancedVolcano(tab_1_all, 
                       x="log2FoldChange", 
                       y="padj", 
                       lab = tab_1_all$Species, 
                       pCutoff = 0.05, 
                       FCcutoff = 0.6,
                       pointSize = 4.0,
                       labSize = 3.0,
                       # boxedLabels = TRUE,
                       title = 'HMO vs Lactose') 



p_1

