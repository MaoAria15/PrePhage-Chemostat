setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################Load sequencing data#######################################

##Import files



source("4_simone.R")


ps <- PSB
alpha_met <- c("Observed") 


standf = function(x, t=total) round(t * (x / sum(x)))
total = median(sample_sums(ps))
PSV.R = transform_sample_counts(ps, standf)

rich <- estimate_richness(PSV.R)
tab <- subset(rich, select = alpha_met)
index <- match(rownames(rich), rownames(sample_data(PSV.R)))
tab$Group <-sample_data(PSV.R)$Group[index]
tab$NEC <-sample_data(PSV.R)$NEC_NO2[index]
tab$Litter <-sample_data(PSV.R)$Litter[index]
#reshape and draw boxplot with ggplot2

tab$Group <- factor(tab$Group,levels = S_groups)
tab$Observed <- round(tab$Observed,2)  ##only for Observed diversity
tab$NEC <- factor(tab$NEC, levels = c("NO","YES"))  

sheet <-list()
###########1. Litter 1-3#####################

tab_sub <- subset(tab, Litter %in% c("1","2","3"))

n1 <-  count(tab_sub$Group == "FVT")[[2,2]]
n2 <-  count(tab_sub$Group == "CVT")[[2,2]]
n3 <-  count(tab_sub$Group == "CVT_MO")[[2,2]]
n4 <-  count(tab_sub$Group == "CON")[[2,2]]

stat_bac<- tab_sub %>%
  wilcox_test(Observed~Group,
              p.adjust.method = "fdr",
              paired = FALSE, 
              alternative = "two.sided",
              detailed = TRUE) 
stat_bac
sheet[["Litter_123"]] <- stat_bac



L_123 <- ggplot(tab_sub, aes(x= Group, y= Observed)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.2) + 
  geom_boxplot(outlier.shape = NA, aes(fill=Group)) +
  geom_point(position = position_jitter(w=0.1, h=0),size=2, aes(colour=NEC)) +
  scale_color_manual(values = NEC_col) +
  labs(color=" NEC")+
  scale_fill_manual(values = simone_cols) +
  scale_x_discrete(labels=c(glue("FVT\nN={n1}"),
                            glue("CVT\nN={n2}"),
                            glue("CVT_MO\nN={n3}"),
                            glue("CON\nN={n4}"))
  )+
  labs(x="", y="Observed alpha diversity", title="16S rRNA") +
  # labs(x="", y="Observed alpha diversity", title="16S rRNA") +
  # stat_pvalue_manual(tat_bac,label = "p.signif", tip.length = 0, size = 6)+ 
  theme_classic() +
  mytheme_beta
L_123

###########2. Litter 45#####################

tab_sub <- subset(tab, Litter %in% c("4","5"))

n1 <-  count(tab_sub$Group == "CVT_MO")[[2,2]]
n2 <-  count(tab_sub$Group == "CON")[[2,2]]

stat_bac<- tab_sub %>%
  wilcox_test(Observed~Group,
              p.adjust.method = "fdr",
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE)
stat_bac

sheet[["Litter_45"]] <- stat_bac

L_45 <- ggplot(tab_sub, aes(x= Group, y= Observed)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.2) + 
  geom_boxplot(outlier.shape = NA, aes(fill=Group)) +
  geom_point(position = position_jitter(w=0.1, h=0),size=2, aes(colour=NEC)) +
  scale_color_manual(values = NEC_col) +
  labs(color=" NEC")+
  scale_fill_manual(values = simone_cols[3:4]) +
  scale_x_discrete(labels=c(glue("CVT_MO\nN={n1}"),
                            glue("CON\nN={n2}"))
  )+
  labs(x="", y="Observed alpha diversity", title="16S rRNA") +
  # labs(x="", y="Observed alpha diversity", title="16S rRNA") +
  # stat_pvalue_manual(tat_bac,label = "p.signif", tip.length = 0, size = 6)+ 
  theme_classic() +
  mytheme_beta
L_45



###########3. Litter 1-3 NEC#####################

tab_sub <- subset(tab, Litter %in% c("1","2","3"))

n1 <-  count(tab_sub$NEC == "NO")[[2,2]]
n2 <-  count(tab_sub$NEC == "YES")[[2,2]]

stat_bac<- tab_sub %>%
  wilcox_test(Observed~NEC,
              p.adjust.method = "fdr",
              paired = FALSE, 
              alternative = "two.sided",
              detailed = TRUE) 
stat_bac
sheet[["Litter_123_nec"]] <- stat_bac

L_123_nec <- ggplot(tab_sub, aes(x= NEC, y= Observed)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.2) + 
  geom_boxplot(outlier.shape = NA, aes(fill=NEC)) +
  geom_point(position = position_jitter(w=0.1, h=0),size=2) +
  scale_fill_manual(values = NEC_col) +
  scale_x_discrete(labels=c(glue("NO\nN={n1}"),
                            glue("YES\nN={n2}"))
  )+
  labs(x="", y="Observed alpha diversity", title="16S rRNA") +
  # labs(x="", y="Observed alpha diversity", title="16S rRNA") +
  stat_pvalue_manual(stat_bac,label = "p", tip.length = 0, size = 3,y.position = c(330))+
  theme_classic() +
  mytheme_beta
L_123_nec

###########4. Litter 45 NEC#####################

tab_sub <- subset(tab, Litter %in% c("4","5"))

n1 <-  count(tab_sub$NEC == "NO")[[2,2]]
n2 <-  count(tab_sub$NEC == "YES")[[2,2]]

# stat_bac<- tab_sub %>%
#   wilcox_test(Observed~NEC,
#               p.adjust.method = "fdr",
#               paired = FALSE, 
#               alternative = "two.sided",
#               detailed = TRUE) 
# stat_bac
# sheet[["Litter_45_nec"]] <- stat_bac

L_45_nec <- ggplot(tab_sub, aes(x= NEC, y= Observed)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.2) + 
  geom_boxplot(outlier.shape = NA, aes(fill=NEC)) +
  geom_point(position = position_jitter(w=0.1, h=0),size=2) +
  scale_fill_manual(values = NEC_col) +
  scale_x_discrete(labels=c(glue("NO\nN={n1}"),
                            glue("YES\nN={n2}"))
  )+
  labs(x="", y="Observed alpha diversity", title="16S rRNA") +
  # labs(x="", y="Observed alpha diversity", title="16S rRNA") +
  # stat_pvalue_manual(tat_bac,label = "p.signif", tip.length = 0, size = 6)+ 
  theme_classic() +
  mytheme_beta
L_45_nec

###########5. Litter 12345 NEC#####################

tab_sub <- subset(tab, Litter %in% c("1","2","3","4","5"))

n1 <-  count(tab_sub$NEC == "NO")[[2,2]]
n2 <-  count(tab_sub$NEC == "YES")[[2,2]]

stat_bac<- tab_sub %>%
  wilcox_test(Observed~NEC,
              p.adjust.method = "fdr",
              paired = FALSE, 
              alternative = "two.sided",
              detailed = TRUE) 
stat_bac
sheet[["Litter_12345_nec"]] <- stat_bac

L_12345_nec <- ggplot(tab_sub, aes(x= NEC, y= Observed)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.2) + 
  geom_boxplot(outlier.shape = NA, aes(fill=NEC)) +
  geom_point(position = position_jitter(w=0.1, h=0),size=2) +
  scale_fill_manual(values = NEC_col) +
  scale_x_discrete(labels=c(glue("NO\nN={n1}"),
                            glue("YES\nN={n2}"))
  )+
  labs(x="", y="Observed alpha diversity", title="16S rRNA") +
  stat_pvalue_manual(stat_bac,label = "p", tip.length = 0, size = 3,y.position = c(380))+
  theme_classic() +
  mytheme_beta
L_12345_nec
write_xlsx(sheet, "Bacteriome/stat_result/Simone/alpha_diversity_Observed.xlsx")



