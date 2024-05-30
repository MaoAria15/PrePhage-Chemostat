setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################Load sequencing data#######################################

##Import files



source("4_simone.R")


ps <- subset_samples(PSV.no.realm, Group %in% S_groups)
alpha_met <- c("Observed","Shannon") 


standf = function(x, t=total) round(t * (x / sum(x)))
total = median(sample_sums(ps))
PSV.R = transform_sample_counts(ps, standf)

rich <- estimate_richness(PSV.R, measures = alpha_met)
tab <- subset(rich, select = alpha_met) 
index <- match(rownames(rich), rownames(sample_data(PSV.R)))
tab$Group <-sample_data(PSV.R)$Group[index]
tab$NEC <-sample_data(PSV.R)$NEC_01[index]
#reshape and draw boxplot with ggplot2

tab$Group <- factor(tab$Group,levels = S_groups)
tab$Shannon <- round(tab$Shannon,2)  ##only for shannon diversity
tab$NEC <- factor(tab$NEC, levels = c("NO","YES"))  

sheet <-list()

################1. Observe by group#################################################

n1 <-  count(tab$Group == "FVT")[[2,2]]
n2 <-  count(tab$Group == "CVT")[[2,2]]
n3 <-  count(tab$Group == "CVT_MO")[[2,2]]
n4 <-  count(tab$Group == "CON")[[2,2]]


stat_bac<- tab %>%
  wilcox_test(Observed~Group,
              p.adjust.method = "fdr",
              paired = FALSE, 
              alternative = "two.sided",
              detailed = TRUE) 
stat_bac
sheet[["Observe_group"]] <- stat_bac

p_alpha <- ggplot(tab, aes(x= Group, y= Observed)) +
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
  labs(x="", y="Observed alpha diversity (Group)", title="16S rRNA") +
  # labs(x="", y="Shannon diversity index", title="16S rRNA") +
  stat_pvalue_manual(stat_bac,label = "p", tip.length = 0, size = 3,
                     y.position = c(90,NA,NA,NA,NA,NA))+
  theme_classic() +
  mytheme_beta
p_alpha

################1. Observe by nec#################################################

n1 <-  count(tab$NEC == "NO")[[2,2]]
n2 <-  count(tab$NEC == "YES")[[2,2]]



stat_bac<- tab %>%
  wilcox_test(Observed~NEC,
              p.adjust.method = "fdr",
              paired = FALSE, 
              alternative = "two.sided",
              detailed = TRUE) 
stat_bac
sheet[["Observe_nec"]] <- stat_bac

p_alpha <- ggplot(tab, aes(x= NEC, y= Observed)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.2) + 
  geom_boxplot(outlier.shape = NA, aes(fill=NEC)) +
  geom_point(position = position_jitter(w=0.1, h=0),size=2) +
  labs(color=" NEC")+
  scale_fill_manual(values = NEC_col) +
  scale_x_discrete(labels=c(glue("NO\nN={n1}"),
                            glue("YES\nN={n2}"))
  )+
  labs(x="", y="Observed alpha diversity (NEC)", title="16S rRNA") +
  # labs(x="", y="Shannon diversity index", title="16S rRNA") +
  # stat_pvalue_manual(stat_bac,label = "p.adj.signif", tip.length = 0, size = 6,
  #                    y.position = c(18,NA,NA,19,NA,NA))+
  theme_classic() +
  mytheme_beta
p_alpha

################3. Shannon by group#################################################

n1 <-  count(tab$Group == "FVT")[[2,2]]
n2 <-  count(tab$Group == "CVT")[[2,2]]
n3 <-  count(tab$Group == "CVT_MO")[[2,2]]
n4 <-  count(tab$Group == "CON")[[2,2]]


stat_bac<- tab %>%
  wilcox_test(Shannon~Group,
              p.adjust.method = "fdr",
              paired = FALSE, 
              alternative = "two.sided",
              detailed = TRUE) 
stat_bac
sheet[["shannon_group"]] <- stat_bac

p_alpha <- ggplot(tab, aes(x= Group, y= Shannon)) +
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
  labs(x="", y="Shannon diversity index (Group)", title="16S rRNA") +
  stat_pvalue_manual(stat_bac,label = "p", tip.length = 0, size = 3,
                     y.position = c(2.5,NA,NA,NA,NA,NA))+
  theme_classic() +
  mytheme_beta
p_alpha

################1. Observe by nec#################################################

n1 <-  count(tab$NEC == "NO")[[2,2]]
n2 <-  count(tab$NEC == "YES")[[2,2]]



stat_bac<- tab %>%
  wilcox_test(Shannon~NEC,
              p.adjust.method = "fdr",
              paired = FALSE, 
              alternative = "two.sided",
              detailed = TRUE) 
stat_bac
sheet[["Shannon_nec"]] <- stat_bac

p_alpha <- ggplot(tab, aes(x= NEC, y= Shannon)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.2) + 
  geom_boxplot(outlier.shape = NA, aes(fill=NEC)) +
  geom_point(position = position_jitter(w=0.1, h=0),size=2) +
  labs(color=" NEC")+
  scale_fill_manual(values = NEC_col) +
  scale_x_discrete(labels=c(glue("NO\nN={n1}"),
                            glue("YES\nN={n2}"))
  )+
  labs(x="", y="Shannon diversity index (NEC)", title="16S rRNA") +
  # labs(x="", y="Shannon diversity index", title="16S rRNA") +
  # stat_pvalue_manual(stat_bac,label = "p.adj.signif", tip.length = 0, size = 6,
  #                    y.position = c(18,NA,NA,19,NA,NA))+
  theme_classic() +
  mytheme_beta
p_alpha
write_xlsx(sheet, "Virome/stat_result/Simone/alpha_diversity.xlsx")
