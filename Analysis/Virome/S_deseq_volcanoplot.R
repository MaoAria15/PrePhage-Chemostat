setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("4_simone.R")
library(EnhancedVolcano)

# level <-"Species"
#############################################Section one : BY GROUP###################################################################################
ps <- subset_samples(PSV.no.realm, Group %in% S_groups)

# ps <- tax_glom(ps, level , NArm = FALSE) #select a level to compare
#draw heatmap 

#Include deseq2 p.adj  taxonomy for rows
#Include Litter Information and Group NEC Information for colors
#load two list of deseq2
tab_1_all <- read.table("Virome/stat_result/Simone/Deseq/FVT_CON.tsv", sep = "\t",header = T, row.names = 1)
tab_2_all <- read.table("Virome/stat_result/Simone/Deseq/CVT_CON.tsv", sep = "\t", header = T, row.names = 1)
tab_3_all <- read.table("Virome/stat_result/Simone/Deseq/CVT_MO_CON.tsv", sep = "\t",header = T, row.names = 1)
tab_4_all <- read.table("Virome/stat_result/Simone/Deseq/CVT_FVT.tsv", sep = "\t", header = T, row.names = 1)
tab_5_all <- read.table("Virome/stat_result/Simone/Deseq/CVT_MO_FVT.tsv", sep = "\t", header = T, row.names = 1)
tab_6_all <- read.table("Virome/stat_result/Simone/Deseq/CVT_CVT_MO.tsv", sep = "\t", header = T, row.names = 1)

tab_1_all$taxa <- paste(rownames(tab_1_all), tab_1_all$Family,sep = "_")
tab_2_all$taxa <- paste(rownames(tab_2_all), tab_2_all$Family,sep = "_")
tab_3_all$taxa <- paste(rownames(tab_3_all), tab_3_all$Family,sep = "_")
tab_4_all$taxa <- paste(rownames(tab_4_all), tab_4_all$Family,sep = "_")
tab_5_all$taxa <- paste(rownames(tab_5_all), tab_5_all$Family,sep = "_")
tab_6_all$taxa <- paste(rownames(tab_6_all), tab_6_all$Family,sep = "_")


###Volcano plot#####################################
#Load packages needed for making volcano plot
#Cut off padj=0.05, log2FoldChange=0.6
p_1<-  EnhancedVolcano(tab_1_all, 
                       x="log2FoldChange", 
                       y="pvalue", 
                       lab = tab_1_all$taxa, 
                       pCutoff = 0.05, 
                       FCcutoff = 0.6,
                       pointSize = 4.0,
                       labSize = 3.0,
                       # boxedLabels = TRUE,
                       title = 'CON vs FVT') 
p_2<-  EnhancedVolcano(tab_2_all, 
                       x="log2FoldChange", 
                       y="pvalue", 
                       lab = tab_2_all$taxa, 
                       pCutoff = 0.05, 
                       FCcutoff = 0.6,
                       pointSize = 4.0,
                       labSize = 3.0,
                       # boxedLabels = TRUE,
                       title = 'CON vs CVT')
p_3<-  EnhancedVolcano(tab_3_all, 
                       x="log2FoldChange", 
                       y="pvalue", 
                       lab = tab_3_all$taxa, 
                       pCutoff = 0.05, 
                       FCcutoff = 0.6,
                       pointSize = 4.0,
                       labSize = 3.0,
                       # boxedLabels = TRUE,
                       title = 'CON vs CVT-MO')
p_4<-  EnhancedVolcano(tab_4_all, 
                       x="log2FoldChange", 
                       y="pvalue", 
                       lab = tab_4_all$taxa, 
                       pCutoff = 0.05, 
                       FCcutoff = 0.6,
                       pointSize = 4.0,
                       labSize = 3.0,
                       # boxedLabels = TRUE,
                       title = 'FVT vs CVT')
p_5<-  EnhancedVolcano(tab_5_all, 
                       x="log2FoldChange", 
                       y="pvalue", 
                       lab =tab_5_all$taxa, 
                       pCutoff = 0.05, 
                       FCcutoff = 0.6,
                       pointSize = 4.0,
                       labSize = 3.0,
                       # boxedLabels = TRUE,
                       title = 'FVT vs CVT-MO')
p_6<-  EnhancedVolcano(tab_6_all, 
                       x="log2FoldChange", 
                       y="pvalue", 
                       lab = tab_6_all$taxa, 
                       pCutoff = 0.05, 
                       FCcutoff = 0.6,
                       pointSize = 4.0,
                       labSize = 3.0,
                       # boxedLabels = TRUE,
                       title = 'CVT-MO vs CVT')

###merge the figure##############################

ggarrange(p_1,p_2,p_3,p_4,p_5,p_6,
          nrow = 2,
          ncol = 3,
          common.legend = T)


#############################################Section one : NEC###################################################################################
ps <- subset_samples(PSV.no.realm, Group %in% S_groups)

# ps <- tax_glom(ps, level , NArm = FALSE) #select a level to compare
#draw heatmap 

#Include deseq2 p.adj  taxonomy for rows
#Include Litter Information and Group NEC Information for colors
#load two list of deseq2
tab_1_all <- read.table("Virome/stat_result/Simone/Deseq/YES_NO_NEC.tsv", sep = "\t",header = T, row.names = 1)

tab_1_all$taxa <- paste(rownames(tab_1_all), tab_1_all$Family,sep = "_")


###Volcano plot#####################################
#Load packages needed for making volcano plot
#Cut off padj=0.05, log2FoldChange=0.6
p_1<-  EnhancedVolcano(tab_1_all, 
                       x="log2FoldChange", 
                       y="pvalue", 
                       lab = tab_1_all$taxa, 
                       pCutoff = 0.05, 
                       FCcutoff = 0.6,
                       pointSize = 4.0,
                       labSize = 3.0,
                       # boxedLabels = TRUE,
                       title = 'No NEC vs NEC') 
p_1
