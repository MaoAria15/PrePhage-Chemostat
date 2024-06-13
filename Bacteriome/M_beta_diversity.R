setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################Load sequencing data#######################################

##Import files



source("4_malene.R")

sheet <- list()
##############Litter 123#############################

ps<-  PSB.CSS
Method="bray"
#############################Weeks compare###########################################


GP.ord <- ordinate(ps, "PCoA", Method)
new_colnames <- gsub("Axis.", "PCoA", colnames(GP.ord$vectors))
colnames(GP.ord$vectors) <- new_colnames
bray.PSB <- phyloseq::distance(ps, method = Method)
# make a data frame from the sample_data
sampledf.PSB <- data.frame(sample_data(ps))
sampledf.PSB <- sampledf.PSB[, c("Group", "Litter")]
# adonis_sub <- adonis2(bray.PSB ~ delivery_mode, method = Method, data = sampledf.PSB, permutations = 999)
adonis_sub <- as.matrix(anova.cca(capscale(bray.PSB ~  sampledf.PSB[,"Group"] + Condition(Litter), sampledf.PSB, dist = Method), permutations = how(nperm=999)))
adonis_R2 <- adonis_sub[1,2]/ (adonis_sub[1,2]+adonis_sub[2,2])
adonis_p <- adonis_sub[1,4]
ps@sam_data$Group <- factor(ps@sam_data$Group, levels = M_groups)
p_beta <- phyloseq::plot_ordination(ps, GP.ord, axes = 1:2,color="Group") + 
  stat_ellipse(geom = "polygon", level = 0.95, fill = NA, linewidth = 1,show.legend = FALSE) +
  geom_point(alpha=2, size=2,aes(color=Group)) +
  scale_color_manual(values= malene_cols)+
  ggtitle(paste( "Bray-curtis\n",
                 "R2=", format(adonis_R2[1],digits = 2),", P=",format(adonis_p[1], digits = 1),"\n")) +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  mytheme
p_beta

#Calculate pairwise p value result - group
metadata <- data.frame(sample_data(ps))
cbn <- combn(x=unique(metadata$Group), m = 2)


p <- c()
r <- c()
f <- c()
for(i in 1:ncol(cbn)){
  ps.subs <- subset_samples(ps, Group %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  permanova_pairwise <- adonis2(phyloseq::distance(ps.subs, method = Method) ~ Group, 
                                data = metadata_sub)
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
  r <- c(r, permanova_pairwise$R2[1])
  f <- c(f, permanova_pairwise$F[1])
}


p_treatment.adj <- round(p.adjust(p, method = "fdr"),digits=3)


t__stat <- data.frame(Group1=c("CON","CON","FVT"),
                            Group2=c("FVT","SPC","SPC"),
                            R2= round(r,digits = 3),
                            F=round(f,digits = 3),
                            p=round(p,digits = 3),
                            p.adj=p_treatment.adj)

t__stat
sheet[["adonis"]] <- t__stat


write_xlsx(sheet,"Bacteriome/stat_result/Malene/beta_diversity.xlsx")
