library(ggplot2)
library(MESS)
library(forcats)
library(tidyverse)
library(psycho)
library(ggpubr)
library(ggsci)
library(rstatix)
library(ampvis2)
library(DESeq2)
library(plyr)
library(cowplot)
library(RVAideMemoire)
library(data.table)
library(microbiome)
library(forestmangr)
library(writexl)
library(viridis)
library(readxl)
library(phyloseq)      # necessary to import the data from Excel file
library(dplyr)        # filter and reformat data frames
library(stringr)
library(vegan)
library(metagenomeSeq)
library(tidyr)
library(RColorBrewer)
library(reshape2)
library(writexl)
library(directlabels)
library(glue)
library(ggpmisc)
library(ComplexHeatmap)
library(magick)
library(colorRamp2)
library(circlize)
library(microDecon)
library(psych)
library(pheatmap)
library(igraph)
library(ggraph)
library(influential)
library(showtext)
library(gridExtra)
library(data.table)
library(EnhancedVolcano)
# library(MicrobiotaProcess)

#Chemostat_pig
simone_cols <- c( "#233462" ,"#0064C7" ,"#16C7DB","#B5743B")
malene_cols <- c("#233462","#B5743B","#16C785")

simone_cols_inoculum <- c( "#233462" ,"#0064C7" ,"#16C7DB","#B5743B","#1b9e77","#d95f02","#7570b3","#e7298a")

S_groups <- c("FVT", "CVT", "CVT_MO","CON")
M_groups <- c("FVT", "CON", "SPC")

NEC_col <- c("black","red")

#Chemostat
HMO_propagation <- c("PV1","PV2","PV4","PV5")
Lactose_propagation <- c("PV3","PV6")
cols_batch <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","pink")
col_donor <- c("#1b9e77","#d95f02","#7570b3","#e7298a")
group_point<- c("Inoculum", "Batch1", "Batch2","Batch3","Chemostat1", "Chemostat2", "Chemostat3", "Chemostat4")
cols_point <- c("darkgrey", "#e31a1c", "#bebada","#fb8072","#80b1d3", "#fdb462", "#b3de69", "#fccde5")
group_subsrate <- c("HMO","Lactose")
col_substrate <- c("red","blue")
col_inoculum <- c("darkgreen")

group_batch <-c("B1","B2","B3","B4","B5","B6")
group_donor <- c("Donor 22 - V.1", "Donor 24 - V.1","Donor 26 - V.1","Donor 27 - V.1")



mytheme_with_x <- theme(text = element_text(size = 8, colour = "Black",family = "Arial"),
                 axis.line=element_line(size=0.5),
                 #panel.border = element_blank(),
                 axis.text=element_text(size = 8, colour = "Black"),
                 axis.ticks=element_line(size=1, colour = "Black"),
                 strip.background = element_rect(colour = "white", fill = "white"),
                 axis.text.x=element_text(size= 8, angle = 0,vjust = 0.6),
                 axis.title = element_text(size = 8, face = "bold"),
                 strip.text.x = element_text(angle = 0, size=8, face = "bold"),
                 legend.text = element_text(size=8),
                 legend.key.size = unit(8, "pt"),
                 legend.title = element_text(size = 8,face = "bold"),
                 title = element_text(size =8, face = "bold")
)
mytheme <- theme(text = element_text(size = 8, colour = "Black",family = "Arial"),
                 axis.line=element_line(size=0.5),
                 #panel.border = element_blank(),
                 axis.text=element_text(size = 8, colour = "Black"),
                 axis.ticks=element_line(size=1, colour = "Black"),
                 axis.ticks.x = element_blank(),
                 strip.background = element_rect(colour = "white", fill = "white"),
                 axis.text.x=element_blank(),
                 axis.title = element_text(size = 8, face = "bold"),
                 strip.text.x = element_text(angle = 30, size=8, face = "bold"),
                 legend.text = element_text(size=8),
                 legend.key.size = unit(8, "pt"),
                 legend.title = element_text(size = 8,face = "bold"),
                 title = element_text(size =8, face = "bold")
)

##Alpha diversity
mytheme_alpha <- theme(text = element_text(size = 8, colour = "Black",family = "Arial"),
                           axis.line=element_line(size=0.5),
                           #panel.border = element_blank(),
                           axis.text.y =element_text(size = 8, colour = "Black"),
                           axis.text.x =element_text(angle = 45,hjust=1,size = 8, colour = "Black"),
                           axis.ticks=element_line(size=1, colour = "Black"),
                           strip.background = element_rect(colour = "white", fill = "white"),
                           axis.title = element_text(size = 8, face = "bold"),
                           strip.text.x = element_text(angle = 30, size=8, face = "bold"),
                           legend.text = element_text(size=8),
                           legend.key.size = unit(8, "pt"),
                           legend.title = element_text(size = 8,face = "bold"),
                           title = element_text(size =8, face = "bold")
)

mytheme_alpha_noy <- theme(text = element_text(size = 8, colour = "Black",family = "Arial"),
                       axis.line=element_line(size=0.5),
                       #panel.border = element_blank(),
                       axis.text=element_text(size = 8, colour = "Black"),
                       axis.ticks=element_line(size=1, colour = "Black"),
                       axis.text.y=element_blank(), 
                       axis.ticks.y=element_blank(), 
                       axis.title.y=element_blank(),
                       axis.line.y=element_blank(),
                       strip.background = element_rect(colour = "white", fill = "white"),
                       # axis.text.x=element_blank(),
                       axis.title = element_text(size = 8, face = "bold"),
                       strip.text.x = element_text(angle = 30, size=8, face = "bold"),
                       legend.text = element_text(size=8),
                       legend.key.size = unit(8, "pt"),
                       legend.title = element_text(size = 8,face = "bold"),
                       title = element_text(size =8, face = "bold")
)
mytheme_beta <- theme(text = element_text(size = 8, colour = "Black",family = "Arial"),
                       axis.line=element_line(size=0.5),
                       #panel.border = element_blank(),
                       axis.text=element_text(size = 8, colour = "Black"),
                       axis.ticks=element_line(size=1, colour = "Black"),
                       strip.background = element_rect(colour = "white", fill = "white"),
                       axis.title = element_text(size = 8, face = "bold"),
                       strip.text.x = element_text( size=8, face = "bold"),
                       legend.text = element_text(size=8),
                       legend.key.size = unit(8, "pt"),
                       legend.title = element_text(size = 8,face = "bold"),
                       title = element_text(size =8, face = "bold")
)
mytheme_abundance <- theme(text = element_text(size = 8, colour = "Black",family = "Arial"),
                      axis.line=element_line(size=0.5),
                      #panel.border = element_blank(),
                      axis.text=element_text(size = 8, colour = "Black"),
                      axis.text.x=element_text(size = 8,colour = "Black"),
                      axis.text.y=element_text(size = 8,colour = "Black"),
                      axis.ticks=element_line(size=1, colour = "Black"),
                      strip.background = element_rect(colour = "white", fill = "white"),
                      axis.title.x  = element_text(size = 8, face = "bold"),
                      strip.text.x = element_text(angle = 0, size=8, face = "bold"),
                      legend.text = element_text(size=8,face = "italic"),
                      legend.key.size = unit(8, "pt"),
                      legend.title = element_text(size = 8,face = "bold"),
                      title = element_text(size =8, face = "bold")
)

mytheme_abundance_noy <- theme(text = element_text(size = 8, colour = "Black",family = "Arial"),
                           axis.line=element_line(size=0.5),
                           axis.text=element_text(size = 8, colour = "Black"),
                           axis.text.x=element_text(size = 8,colour = "Black"),
                           axis.text.y=element_blank(), 
                           axis.ticks.y=element_blank(), 
                           axis.title.y=element_blank(),
                           axis.line.y=element_blank(),
                           axis.ticks=element_line(size=1, colour = "Black"),
                           strip.background = element_rect(colour = "white", fill = "white"),
                           axis.title = element_text(size = 8, face = "bold"),
                           strip.text.x = element_text(angle = 0, size=8, face = "bold"),
                           legend.text = element_text(size=8,face = "italic"),
                           legend.key.size = unit(8, "pt"),
                           legend.title = element_text(size = 8,face = "bold"),
                           title = element_text(size =8, face = "bold"))


mytheme_abundance_yvertival <- theme(text = element_text(size = 8, colour = "Black",family = "Arial"),
                           axis.line=element_line(size=0.5),
                           #panel.border = element_blank(),
                           axis.text=element_text(size = 8, colour = "Black"),
                           axis.text.x=element_text(size = 8,angle=90,colour = "Black"),
                           axis.text.y=element_text(size = 8,colour = "Black"),
                           axis.ticks=element_line(size=1, colour = "Black"),
                           strip.background = element_rect(colour = "white", fill = "white"),
                           axis.title.x  = element_text(size = 8, face = "bold"),
                           strip.text.x = element_text(angle = 0, size=8, face = "bold"),
                           legend.text = element_text(size=8,face = "italic"),
                           legend.key.size = unit(8, "pt"),
                           legend.title = element_text(size = 8,face = "bold"),
                           title = element_text(size =8, face = "bold")
)





mytheme_his <-  theme(text = element_text(size = 8, colour = "Black",family = "Arial"),
                      axis.line=element_line(size=0.5),
                      #panel.border = element_blank(),
                      axis.text=element_text(size = 8, colour = "Black"),
                      axis.ticks=element_line(size=1, colour = "Black"),
                      strip.background = element_rect(colour = "white", fill = "white"),
                      axis.text.x=element_text(size= 8, angle = 0,vjust = 0.6),
                      axis.title = element_text(size = 8, face = "bold"),
                      strip.text.x = element_text(angle = 30, size=8, face = "bold"),
                      legend.text = element_text(size=8),
                      legend.key.size = unit(8, "pt"),
                      legend.title = element_text(size = 8,face = "bold"),
                      title = element_text(size =8, face = "bold"),
                      legend.background = element_rect(color = "black",linewidth = 1))


filter_and_replace <- function(data, threshold = 0.05) {
  data %>%
    filter(p < threshold) %>%
    mutate(p_signif = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ ""
    )) %>%
    select(-p)  # Remove the original p-value column
}


