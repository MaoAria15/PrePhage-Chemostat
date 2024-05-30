#Set working directory to script directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("4_simone.R")

set.seed(19950915)
ps<- PSV.no.Realm.CSS

#Calculate bray curtis dissimilarity within each group
# Divide phyloseq into each group
FVT <- subset_samples(ps,Group == "FVT")
CVT <- subset_samples(ps,Group == "CVT")
CVT_MO <- subset_samples(ps,Group == "CVT_MO")
CON <- subset_samples(ps,Group == "CON")

# calculate distance
dis_FVT <- as.matrix(phyloseq::distance(FVT, method = "bray"))
dis_CVT <- as.matrix(phyloseq::distance(CVT, method = "bray"))
dis_CVT_MO <- as.matrix(phyloseq::distance(CVT_MO, method = "bray"))
dis_CON <- as.matrix(phyloseq::distance(CON, method = "bray"))

sub_design_FVT <- data.frame(sample_data(FVT))
sub_design_CVT <- data.frame(sample_data(CVT))
sub_design_CVT_MO <- data.frame(sample_data(CVT_MO))
sub_design_CON <- data.frame(sample_data(CON))

dis_FVT <- dis_FVT[rownames(sub_design_FVT),rownames(sub_design_FVT)]
dis_CVT <- dis_CVT[rownames(sub_design_CVT),rownames(sub_design_CVT)]
dis_CVT_MO <- dis_CVT_MO[rownames(sub_design_CVT_MO),rownames(sub_design_CVT_MO)]
dis_CON <- dis_CON[rownames(sub_design_CON),rownames(sub_design_CON)]

# Get the upper triangle indices
# FVT
diag(dis_FVT) <- NA
dis_FVT[upper.tri(dis_FVT)] <-NA
long_FVT <- melt(dis_FVT, na.rm = TRUE)
FVT_var<-long_FVT$value

#CVT
diag(dis_CVT) <- NA
dis_CVT[upper.tri(dis_CVT)] <-NA
long_CVT <- melt(dis_CVT, na.rm = TRUE)
CVT_var <- long_CVT$value
#dis_CVT_MO
diag(dis_CVT_MO) <- NA
dis_CVT_MO[upper.tri(dis_CVT_MO)] <-NA
long_CVT_MO <- melt(dis_CVT_MO, na.rm = TRUE)
CVT_MO_var<-long_CVT_MO$value
#CON
diag(dis_CON) <- NA
dis_CON[upper.tri(dis_CON)] <-NA
long_CON <- melt(dis_CON, na.rm = TRUE)
CON_var<-long_CON$value


#Combine dissimilarity into long datafram


combined_dis <- data.frame(group = rep(c("FVT", "CVT", "CVT_MO", "CON"), times = c(length(FVT_var), length(CVT_var), length(CVT_MO_var),length(CON_var))),
                           dissimilarity = c(FVT_var, CVT_var, CVT_MO_var,CON_var))



combined_dis$group <- factor(combined_dis$group ,levels = c("FVT", "CVT", "CVT_MO", "CON"))




stat <- wilcox_test( dissimilarity ~ group, 
                     data = combined_dis,
                     p.adjust.method = "fdr" )

stat

ggplot(combined_dis, aes(x = group, y = dissimilarity))+
  stat_boxplot(geom ='errorbar', linetype=1, width=0.2) + 
  geom_boxplot(outlier.shape = NA, aes(fill=group)) +
  scale_fill_manual(values = simone_cols) +
  labs(x="", y="Bray-Curtis dissimilarity") +
  # stat_pvalue_manual(stat2,label = "p", tip.length = 0, size = 4,
  #                    y.position = c(NA,NA,0.55,NA,0.6,0.65))+
  theme_classic() +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1))+
  theme(text = element_text(size = 10, colour = "Black"),
        axis.line=element_line(size=0.5),
        #panel.border = element_blank(),
        axis.text=element_text(size = 10, colour = "Black"),
        axis.ticks=element_line(size=1, colour = "Black"),
        strip.background = element_rect(colour = "white", fill = "white"),
        axis.text.x=element_text(size= 10, angle = 0),
        axis.title = element_text(size = 12, face = "bold"),
        strip.text.x = element_text(angle = 0, size=12, face = "bold"),
        legend.text = element_text(size=10),
        legend.key.size = unit(10, "pt"),
        legend.title = element_text(size = 12,face = "bold"),
        title = element_text(size =14, face = "bold")
  ) 
write.csv(stat,file = "Virome/stat_result/Simone/heterogenity_bray.csv")


