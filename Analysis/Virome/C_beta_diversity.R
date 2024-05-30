setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################Load sequencing data#######################################

##Import files



# source("4_Chemostat.R")
source("5_Chemostat_deduct_media.R") #deduct media
sheet <- list()
##############HMO#############################

ps<-  subset_samples(PSV.no.Realm.CSS,Substrate %in% "HMO")
ps<-  subset_samples(ps,Sample_time_point %in% group_point[2:8])
Method="bray"



GP.ord <- ordinate(ps, "PCoA", Method)
new_colnames <- gsub("Axis.", "PCoA", colnames(GP.ord$vectors))
colnames(GP.ord$vectors) <- new_colnames
bray.PSB <- phyloseq::distance(ps, method = Method)
# make a data frame from the sample_data
sampledf.PSB <- data.frame(sample_data(ps))
sampledf.PSB <- sampledf.PSB[, c("Sample_time_point", "Batch")]
# adonis_sub <- adonis2(bray.PSB ~ delivery_mode, method = Method, data = sampledf.PSB, permutations = 999)
adonis_sub <- as.matrix(anova.cca(capscale(bray.PSB ~  sampledf.PSB[,"Sample_time_point"] + Condition(Batch), sampledf.PSB, dist = Method), permutations = how(nperm=999)))
adonis_R2 <- adonis_sub[1,2]/ (adonis_sub[1,2]+adonis_sub[2,2])
adonis_p <- adonis_sub[1,4]
ps@sam_data$Sample_time_point <- factor(ps@sam_data$Sample_time_point, levels = group_point[2:8])
p_HMO <- phyloseq::plot_ordination(ps, GP.ord, axes = 1:2,color="Sample_time_point") + 
  stat_ellipse(geom = "polygon", level = 0.95, fill = NA, linewidth = 1,show.legend = FALSE) +
  geom_point(alpha=2, size=2,aes(color=Sample_time_point)) +
  scale_color_manual(values= cols_point[2:8])+
  ggtitle(paste( "HMO\n",
                 "R2=", format(adonis_R2[1],digits = 2),", P=",format(adonis_p[1], digits = 1),"\n")) +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  mytheme
p_HMO

#Calculate pairwise p value result - group
metadata <- data.frame(sample_data(ps))
cbn <- combn(x=unique(metadata$Sample_time_point), m = 2)


p <- c()
r <- c()
f <- c()
for(i in 1:ncol(cbn)){
  ps.subs <- subset_samples(ps, Sample_time_point %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  permanova_pairwise <- adonis2(phyloseq::distance(ps.subs, method = Method) ~ Sample_time_point, 
                                data = metadata_sub)
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
  r <- c(r, permanova_pairwise$R2[1])
  f <- c(f, permanova_pairwise$F[1])
}


p_treatment.adj <- round(p.adjust(p, method = "fdr"),digits=3)


t__HMO <- data.frame(Group1=c("Batch1","Batch1","Batch1","Batch1","Batch1","Batch1","Batch2","Batch2","Batch2","Batch2","Batch2","Batch3","Batch3","Batch3","Batch3","Chemostat1","Chemostat1","Chemostat1","Chemostat2","Chemostat2","Chemostat3"),
                     Group2=c("Batch2","Batch3","Chemostat1","Chemostat2","Chemostat3","Chemostat4","Batch3","Chemostat1","Chemostat2","Chemostat3","Chemostat4","Chemostat1","Chemostat2","Chemostat3","Chemostat4","Chemostat2","Chemostat3","Chemostat4","Chemostat3","Chemostat4","Chemostat4"),
                     R2= round(r,digits = 3),
                     F=round(f,digits = 3),
                     p=round(p,digits = 3),
                     p.adj=p_treatment.adj)

t__HMO
sheet[["HMO"]] <- t__HMO
##############Lactose#############################


ps<-  subset_samples(PSV.no.Realm.CSS,Substrate %in% "Lactose")
ps<-  subset_samples(ps,Sample_time_point %in% group_point[2:8])
Method="bray"



GP.ord <- ordinate(ps, "PCoA", Method)
new_colnames <- gsub("Axis.", "PCoA", colnames(GP.ord$vectors))
colnames(GP.ord$vectors) <- new_colnames
bray.PSB <- phyloseq::distance(ps, method = Method)
# make a data frame from the sample_data
sampledf.PSB <- data.frame(sample_data(ps))
sampledf.PSB <- sampledf.PSB[, c("Sample_time_point", "Batch")]
# adonis_sub <- adonis2(bray.PSB ~ delivery_mode, method = Method, data = sampledf.PSB, permutations = 999)
adonis_sub <- as.matrix(anova.cca(capscale(bray.PSB ~  sampledf.PSB[,"Sample_time_point"] + Condition(Batch), sampledf.PSB, dist = Method), permutations = how(nperm=999)))
adonis_R2 <- adonis_sub[1,2]/ (adonis_sub[1,2]+adonis_sub[2,2])
adonis_p <- adonis_sub[1,4]
ps@sam_data$Sample_time_point <- factor(ps@sam_data$Sample_time_point, levels = group_point[2:8])
p_Lactose <- phyloseq::plot_ordination(ps, GP.ord, axes = 1:2,color="Sample_time_point") + 
  stat_ellipse(geom = "polygon", level = 0.95, fill = NA, linewidth = 1,show.legend = FALSE) +
  geom_point(alpha=2, size=2,aes(color=Sample_time_point)) +
  scale_color_manual(values= cols_point[2:8])+
  ggtitle(paste( "Lactose\n",
                 "R2=", format(adonis_R2[1],digits = 2),", P=",format(adonis_p[1], digits = 1),"\n")) +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  mytheme
p_Lactose




# Calculate pairwise p value result - group
#Note because of low sample size the result maybe not be reliable
metadata <- data.frame(sample_data(ps))
cbn <- combn(x=unique(metadata$Sample_time_point), m = 2)


p <- c()
r <- c()
f <- c()
for(i in 1:ncol(cbn)){
  ps.subs <- subset_samples(ps, Sample_time_point %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  permanova_pairwise <- adonis2(phyloseq::distance(ps.subs, method = Method) ~ Sample_time_point,
                                data = metadata_sub,permutations = how(nperm=999))
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
  r <- c(r, permanova_pairwise$R2[1])
  f <- c(f, permanova_pairwise$F[1])
}


p_treatment.adj <- round(p.adjust(p, method = "fdr"),digits=3)


t__Lactose <- data.frame(Group1=c("Batch1","Batch1","Batch1","Batch1","Batch1","Batch1","Batch2","Batch2","Batch2","Batch2","Batch2","Batch3","Batch3","Batch3","Batch3","Chemostat1","Chemostat1","Chemostat1","Chemostat2","Chemostat2","Chemostat3"),
                         Group2=c("Batch2","Batch3","Chemostat1","Chemostat2","Chemostat3","Chemostat4","Batch3","Chemostat1","Chemostat2","Chemostat3","Chemostat4","Chemostat1","Chemostat2","Chemostat3","Chemostat4","Chemostat2","Chemostat3","Chemostat4","Chemostat3","Chemostat4","Chemostat4"),
                         R2= round(r,digits = 3),
                         F=round(f,digits = 3),
                         p=round(p,digits = 3),
                         p.adj=p_treatment.adj)

t__Lactose
sheet[["Lactose"]] <- t__Lactose
##############All samples#############################

ps<-  subset_samples(PSV.no.Realm.CSS,Sample_time_point %in% group_point)
Method="bray"



GP.ord <- ordinate(ps, "PCoA", Method)
new_colnames <- gsub("Axis.", "PCoA", colnames(GP.ord$vectors))
colnames(GP.ord$vectors) <- new_colnames
bray.PSB <- phyloseq::distance(ps, method = Method)
# make a data frame from the sample_data
sampledf.PSB <- data.frame(sample_data(ps))
sampledf.PSB <- sampledf.PSB[, c("Sample_time_point", "Batch")]
# adonis_sub <- adonis2(bray.PSB ~ delivery_mode, method = Method, data = sampledf.PSB, permutations = 999)
adonis_sub <- as.matrix(anova.cca(capscale(bray.PSB ~  sampledf.PSB[,"Sample_time_point"] + Condition(Batch), sampledf.PSB, dist = Method), permutations = how(nperm=999)))
adonis_R2 <- adonis_sub[1,2]/ (adonis_sub[1,2]+adonis_sub[2,2])
adonis_p <- adonis_sub[1,4]
ps@sam_data$Sample_time_point <- factor(ps@sam_data$Sample_time_point, levels = group_point)
ps@sam_data$Substrate <- factor(ps@sam_data$Substrate, levels = c("Inoculum","HMO","Lactose"))
p_all_sample <- phyloseq::plot_ordination(ps, GP.ord, axes = 1:2,color="Sample_time_point") + 
  stat_ellipse(geom = "polygon", level = 0.95, fill = NA, linewidth = 1,show.legend = FALSE) +
  geom_point(alpha=2, size=2,aes(color=Sample_time_point)) +
  scale_color_manual(values= cols_point)+
  facet_wrap("Substrate")+
  ggtitle(paste( "All sample\n",
                 "R2=", format(adonis_R2[1],digits = 2),", P=",format(adonis_p[1], digits = 1),"\n")) +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  mytheme_beta
p_all_sample




# Calculate pairwise p value result - group
#Note because of low sample size the result maybe not be reliable
metadata <- data.frame(sample_data(ps))
cbn <- combn(x=unique(metadata$Sample_time_point), m = 2)


p <- c()
r <- c()
f <- c()
for(i in 1:ncol(cbn)){
  ps.subs <- subset_samples(ps, Sample_time_point %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  permanova_pairwise <- adonis2(phyloseq::distance(ps.subs, method = Method) ~ Sample_time_point,
                                data = metadata_sub,permutations = how(nperm=999))
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
  r <- c(r, permanova_pairwise$R2[1])
  f <- c(f, permanova_pairwise$F[1])
}


p_treatment.adj <- round(p.adjust(p, method = "fdr"),digits=3)


t__all_sample <- data.frame(Group1=c("Inoculum","Inoculum","Inoculum","Inoculum","Inoculum","Inoculum","Inoculum","Batch1","Batch1","Batch1","Batch1","Batch1","Batch1","Batch2","Batch2","Batch2","Batch2","Batch2","Batch3","Batch3","Batch3","Batch3","Chemostat1","Chemostat1","Chemostat1","Chemostat2","Chemostat2","Chemostat3"),
                            Group2=c("Batch1","Batch2","Batch3","Chemostat1","Chemostat2","Chemostat3","Chemostat4","Batch2","Batch3","Chemostat1","Chemostat2","Chemostat3","Chemostat4","Batch3","Chemostat1","Chemostat2","Chemostat3","Chemostat4","Chemostat1","Chemostat2","Chemostat3","Chemostat4","Chemostat2","Chemostat3","Chemostat4","Chemostat3","Chemostat4","Chemostat4"),
                            R2= round(r,digits = 3),
                            F=round(f,digits = 3),
                            p=round(p,digits = 3),
                            p.adj=p_treatment.adj)

t__all_sample
sheet[["all_sample"]] <- t__all_sample


#####################################################################################

p_all_sample

ggarrange(p_HMO,p_Lactose,
          nrow = 1,
          ncol = 2,
          common.legend = T,
          legend = "right") #800*400


###################################################################################
write_xlsx(sheet, "Virome/stat_result/Chemostat/beta_diversity(deduct_media).xlsx")
# write_xlsx(sheet, "Virome/stat_result/Chemostat/beta_diversity.xlsx")
