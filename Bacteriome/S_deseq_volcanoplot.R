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
tab_1_all <- read.table("Bacteriome/stat_result/Simone/Deseq/Litter123_species/FVT_CON.tsv", sep = "\t",header = T, row.names = 1)
tab_2_all <- read.table("Bacteriome/stat_result/Simone/Deseq/Litter123_species/CVT_CON.tsv", sep = "\t", header = T, row.names = 1)
tab_3_all <- read.table("Bacteriome/stat_result/Simone/Deseq/Litter123_species/CVT_MO_CON.tsv", sep = "\t",header = T, row.names = 1)
tab_4_all <- read.table("Bacteriome/stat_result/Simone/Deseq/Litter123_species/CVT_FVT.tsv", sep = "\t", header = T, row.names = 1)
tab_5_all <- read.table("Bacteriome/stat_result/Simone/Deseq/Litter123_species/CVT_MO_FVT.tsv", sep = "\t", header = T, row.names = 1)
tab_6_all <- read.table("Bacteriome/stat_result/Simone/Deseq/Litter123_species/CVT_CVT_MO.tsv", sep = "\t", header = T, row.names = 1)
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
                                  title = 'CON vs CVT')
p_3<-  EnhancedVolcano(tab_3_all, 
                                  x="log2FoldChange", 
                                  y="padj", 
                                  lab = tab_3_all$Species, 
                                  pCutoff = 0.05, 
                                  FCcutoff = 0.6,
                                  pointSize = 4.0,
                                  labSize = 3.0,
                                  # boxedLabels = TRUE,
                                  title = 'CON vs CVT-MO')
p_4<-  EnhancedVolcano(tab_4_all, 
                                  x="log2FoldChange", 
                                  y="padj", 
                                  lab = tab_4_all$Species, 
                                  pCutoff = 0.05, 
                                  FCcutoff = 0.6,
                                  pointSize = 4.0,
                                  labSize = 3.0,
                                  # boxedLabels = TRUE,
                                  title = 'FVT vs CVT')
p_5<-  EnhancedVolcano(tab_5_all, 
                                  x="log2FoldChange", 
                                  y="padj", 
                                  lab =tab_5_all$Species, 
                                  pCutoff = 0.05, 
                                  FCcutoff = 0.6,
                                  pointSize = 4.0,
                                  labSize = 3.0,
                                  # boxedLabels = TRUE,
                                  title = 'FVT vs CVT-MO')
p_6<-  EnhancedVolcano(tab_6_all, 
                                  x="log2FoldChange", 
                                  y="padj", 
                                  lab = tab_6_all$Species, 
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



##################Section TWO : litter 45###################################################################################
ps <- subset_samples(PSB, Litter %in% c("4","5"))

ps <- tax_glom(ps, level , NArm = FALSE) #select a level to compare
#draw heatmap 

#Include deseq2 p.adj  taxonomy for rows
#Include Litter Information and Group NEC Information for colors
#load two list of deseq2
tab_1_all <- read.table("Bacteriome/stat_result/Simone/Deseq/Litter45_species/CVT_MO_CON.tsv", sep = "\t",header = T, row.names = 1)
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
                       title = 'CVT-MO vs CON') 
p_1

