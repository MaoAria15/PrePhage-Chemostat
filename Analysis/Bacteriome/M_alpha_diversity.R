setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################Load sequencing data#######################################

##Import files



source("4_malene.R")


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
#reshape and draw boxplot with ggplot2

tab$Group <- factor(tab$Group,levels = M_groups)
tab$Observed <- round(tab$Observed,2)  ##only for shannon diversity
tab$NEC <- factor(tab$NEC, levels = c("NO","YES"))  

sheet <-list()
##draw figures
n1 <-  count(tab$Group == "FVT")[[2,2]]
n2 <-  count(tab$Group == "CON")[[2,2]]
n3 <-  count(tab$Group == "SPC")[[2,2]]


stat_bac<- tab %>%
  wilcox_test(Observed~Group,
              p.adjust.method = "fdr",
              paired = FALSE, 
              alternative = "two.sided",
              detailed = TRUE) 
stat_bac
sheet[["stat_alpha"]] <- stat_bac

p_alpha <- ggplot(tab, aes(x= Group, y= Observed)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.2) + 
  geom_boxplot(outlier.shape = NA, aes(fill=Group)) +
  geom_point(position = position_jitter(w=0.1, h=0),size=2, aes(colour=NEC)) +
  scale_color_manual(values = NEC_col) +
  labs(color=" NEC")+
  scale_fill_manual(values = malene_cols) +
  scale_x_discrete(breaks=M_groups,
                   labels=c(glue("FVT\nN={n1}"),
                            glue("CON\nN={n2}"),
                            glue("SPC\nN={n3}")
                   )
  ) +
  scale_y_continuous(limits = c(170, 350))+
  labs(x="", y="Observed alpha diversity", title="16S rRNA") +
  # labs(x="", y="Shannon diversity index", title="16S rRNA") +
  stat_pvalue_manual(stat_bac,label = "p.adj", tip.length = 0, size = 4,
                     y.position = c(350,NA,345))+
  theme_classic() +
  mytheme_beta
p_alpha



write_xlsx(sheet, "Bacteriome/stat_result/Malene/alpha_diversity.xlsx")
