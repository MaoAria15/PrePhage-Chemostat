setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# source("4_Chemostat.R")
source("5_Chemostat_deduct_media.R") #dedect media

library(EnhancedVolcano)


#############################################Section one : litter 123###################################################################################
ps <- subset_samples(PSV.no.realm, Sample_time_point %in% c("Chemostat4"))
#draw heatmap 

#Include deseq2 p.adj  taxonomy for rows
#load two list of deseq2
# tab_1_all <- read.table("Virome/stat_result/Chemostat/Deseq/Lactose_HMO_chemostat4.tsv", sep = "\t",header = T, row.names = 1)
tab_1_all <- read.table("Virome/stat_result/Chemostat/Deseq/Lactose_HMO_chemostat4(deduct_media).tsv", sep = "\t",header = T, row.names = 1)

tab_1_all$taxa <- paste(rownames(tab_1_all), tab_1_all$Family,sep = "_")
###Volcano plot#####################################
#Load packages needed for making volcano plot
#Cut off padj=0.05, log2FoldChange=0.6
p_1<-  EnhancedVolcano(tab_1_all, 
                       x="log2FoldChange", 
                       y="padj", 
                       lab = tab_1_all$taxa, 
                       pCutoff = 0.05, 
                       FCcutoff = 0.6,
                       pointSize = 4.0,
                       labSize = 3.0,
                       # boxedLabels = TRUE,
                       title = 'HMO vs Lactose') 



p_1

