setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("4_malene.R")
library(EnhancedVolcano)

level <-"Species"


ps <- tax_glom(PSB, level , NArm = FALSE) #select a level to compare
#draw heatmap 

#Include deseq2 p.adj  taxonomy for rows
#Include Litter Information and Group NEC Information for colors
#load two list of deseq2
tab_1_all <- read.table("Bacteriome/stat_result/Malene/Deseq/FVT_CON.tsv", sep = "\t",header = T, row.names = 1)
tab_2_all <- read.table("Bacteriome/stat_result/Malene/Deseq/SPC_CON.tsv", sep = "\t", header = T, row.names = 1)
tab_3_all <- read.table("Bacteriome/stat_result/Malene/Deseq/SPC_FVT.tsv", sep = "\t",header = T, row.names = 1)

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
                       title = 'CON vs FVT') 
p_2<-  EnhancedVolcano(tab_2_all, 
                       x="log2FoldChange", 
                       y="padj", 
                       lab = tab_2_all$Species, 
                       pCutoff = 0.05, 
                       FCcutoff = 0.6,
                       pointSize = 4.0,
                       labSize = 3.0,
                       # boxedLabels = TRUE,
                       title = 'CON vs SPC')
p_3<-  EnhancedVolcano(tab_3_all, 
                       x="log2FoldChange", 
                       y="padj", 
                       lab = tab_3_all$Species, 
                       pCutoff = 0.05, 
                       FCcutoff = 0.6,
                       pointSize = 4.0,
                       labSize = 3.0,
                       # boxedLabels = TRUE,
                       title = 'FVT vs SPC')


###merge the figure##############################

ggarrange(p_1,p_2,p_3,
          nrow = 1,
          ncol = 3,
          common.legend = T)



