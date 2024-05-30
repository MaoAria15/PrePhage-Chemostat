setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("4_simone.R")
library(EnhancedVolcano)

level <-"Species"
#############################################Section one : litter 123###################################################################################
ps <- subset_samples(PSB, Litter %in% c("1","2","3"))

ps <- tax_glom(ps, level , NArm = FALSE) #select a level to compare
#draw heatmap 

#Include deseq2 p.adj  taxonomy for rows
#Include Litter Information and Group NEC Information for colors
#load two list of deseq2
tab_1_all <- read.table("Bacteriome/stat_result/Simone/Deseq/Litter123_species/YES_NO_litter123.tsv", sep = "\t",header = T, row.names = 1)

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
                       title = 'No NEC vs NEC (Litter 123)') 


p_1

##################Section TWO : litter 45###################################################################################
ps <- subset_samples(PSB, Litter %in% c("4","5"))

ps <- tax_glom(ps, level , NArm = FALSE) #select a level to compare
#draw heatmap 

#Include deseq2 p.adj  taxonomy for rows
#Include Litter Information and Group NEC Information for colors
#load two list of deseq2
tab_1_all <- read.table("Bacteriome/stat_result/Simone/Deseq/Litter45_species/YES_NO_litter45.tsv", sep = "\t",header = T, row.names = 1)
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
                       title = 'No NEC vs NEC (Litter 45)') 
p_1

##################Section TWO : litter 45###################################################################################
ps <- subset_samples(PSB, Litter %in% c("1","2","3","4","5"))

ps <- tax_glom(ps, level , NArm = FALSE) #select a level to compare
#draw heatmap 

#Include deseq2 p.adj  taxonomy for rows
#Include Litter Information and Group NEC Information for colors
#load two list of deseq2
tab_1_all <- read.table("Bacteriome/stat_result/Simone/Deseq/YES_NO_litter12345.tsv", sep = "\t",header = T, row.names = 1)
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
                       title = 'No NEC vs NEC (Litter 12345)') 
p_1
