setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################Load sequencing data#######################################

##Import files



source("4_Chemostat.R")
# source("5_Chemostat_deduct_media.R") #Deduct media

ps <- PSV.no.realm
alpha_met <- c("Observed") 


standf = function(x, t=total) round(t * (x / sum(x)))
total = median(sample_sums(ps))
PSV.R = transform_sample_counts(ps, standf)

rich <- estimate_richness(PSV.R)
tab <- subset(rich, select = alpha_met)
index <- match(rownames(rich), rownames(sample_data(PSV.R)))
tab$Substrate <-sample_data(PSV.R)$Substrate[index]
tab$Sample_time_point <-sample_data(PSV.R)$Sample_time_point[index]

#reshape and draw boxplot with ggplot2
tab_long <- melt(tab,id.vars = c("Substrate","Sample_time_point"),
                 variable.name = "Parameter", value.name = "Abundance")
tab_long


tab_long$Substrate <- factor(tab_long$Substrate,levels = group_subsrate)
tab_long$Abundance <- round(tab_long$Abundance,2)  ##only for shannon diversity
#########Section 1: divided by sample time point###########
####################Function setup##########################
count_number <- function(df)
{
  vector<-c()
  vector[1] <-  count(df$Substrate == "HMO")[[2,2]]
  vector[2] <-  count(df$Substrate == "Lactose")[[2,2]]
  return(vector)
}


stat_single_compare <- function(df){
  stat_bac<- df %>%
    wilcox_test(Abundance~Substrate,
                # p.adjust.method = "fdr",
                paired = FALSE, 
                alternative = "two.sided",
                detailed = TRUE) 
  return(stat_bac)
}



ggplot_boxplot_abundance <- function(df)
{
  ggplot(df, aes(x= Substrate, y= Abundance))+ 
    scale_fill_manual(values = col_substrate) +
    stat_boxplot(geom ='errorbar', linetype=1, width=0.2) + 
    geom_boxplot(outlier.shape = NA, aes(fill=Substrate)) +
    geom_jitter(show.legend=FALSE, width=0.25, shape=21, fill="black") +
    scale_x_discrete(breaks=group_subsrate,
                     labels=c(glue("HMO\nN={counts[1]}"),
                              glue("Lactose\nN={counts[2]}"))) +
    # scale_y_continuous(limits = c(10, 45))+
    theme_classic() +
    mytheme_alpha_noy
}

sheet <- list()
######################Plot drawing and statistic cal##############################################################################
# 0. Inoculum####

variable <- "Inoculum"
df<-subset(tab_long, Sample_time_point %in% variable)
counts <- count(df$Sample_time_point == variable)[[1,2]]

p_inoculum  <-  ggplot(df, aes(x= Sample_time_point, y= Abundance))+ 
  stat_boxplot(geom ='errorbar', linetype=1, width=0.2) + 
  geom_boxplot(outlier.shape = NA, aes(fill=Substrate)) +
  scale_fill_manual(values = cols_point[1]) + 
  geom_jitter(show.legend=FALSE, width=0.25, shape=21, fill="black") +
  scale_x_discrete(labels=c(glue("Inoculum\nN={counts[1]}"))) +
  theme_classic() +
  mytheme_beta+
  labs(x="", y="Observe alpha diversity", title=variable)+
  # scale_y_continuous(limits = c(10, 45))+
  theme(plot.title = element_text(hjust = 0.5))

p_inoculum

# 1. Batch1####

variable <- "Batch1"
df<-subset(tab_long, Sample_time_point %in% variable)
stat <- stat_single_compare(df)
stat
counts <- count_number(df)

B_1 <-  ggplot_boxplot_abundance(df)+
  labs(x="",  title=variable)+
  theme(plot.title = element_text(hjust = 0.5))

B_1


# 2. Batch2####

variable <- "Batch2"
df<-subset(tab_long, Sample_time_point %in% variable)
stat <- stat_single_compare(df)
stat
counts <- count_number(df)

B_2 <-  ggplot_boxplot_abundance(df)+
  labs(x="",  title=variable)+
  theme(plot.title = element_text(hjust = 0.5))

B_2



# 3. Batch3####

variable <- "Batch3"
df<-subset(tab_long, Sample_time_point %in% variable)
stat <- stat_single_compare(df)
stat
counts <- count_number(df)

B_3 <-  ggplot_boxplot_abundance(df)+
  labs(x="",  title=variable)+
  theme(plot.title = element_text(hjust = 0.5))

B_3



# 4. Chemostat1####

variable <- "Chemostat1"
df<-subset(tab_long, Sample_time_point %in% variable)
stat <- stat_single_compare(df)
stat
counts <- count_number(df)

C_1 <-  ggplot_boxplot_abundance(df)+
  labs(x="",  title=variable)+
  theme(plot.title = element_text(hjust = 0.5))

C_1



# 5. Chemostat2####

variable <- "Chemostat2"
df<-subset(tab_long, Sample_time_point %in% variable)
stat <- stat_single_compare(df)
stat
counts <- count_number(df)

C_2 <-  ggplot_boxplot_abundance(df)+
  labs(x="",  title=variable)+
  theme(plot.title = element_text(hjust = 0.5))

C_2


# 6. Chemostat3####

variable <- "Chemostat3"
df<-subset(tab_long, Sample_time_point %in% variable)
stat <- stat_single_compare(df)
stat
counts <- count_number(df)

C_3 <-  ggplot_boxplot_abundance(df)+
  labs(x="",  title=variable)+
  theme(plot.title = element_text(hjust = 0.5))

C_3


# 7. Chemostat4####

variable <- "Chemostat4"
df<-subset(tab_long, Sample_time_point %in% variable)
stat <- stat_single_compare(df)
stat
counts <- count_number(df)

C_4 <-  ggplot_boxplot_abundance(df)+
  labs(x="",  title=variable)+
  theme(plot.title = element_text(hjust = 0.5))

C_4




################merge all figures together################

ggarrange(p_inoculum,
          B_1,B_2,B_3,
          C_1,C_2,C_3,C_4,
          nrow = 1,
          ncol = 8,
          common.legend = T,
          widths = c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5)) # For export svg 900x400 size

#########Section 2: divided by substrate###########


####################Function setup##########################
count_number <- function(df)
{
  vector<-c()
  vector[1] <-  count(df$Sample_time_point == "Batch1")[[2,2]]
  vector[2] <-  count(df$Sample_time_point == "Batch2")[[2,2]]
  vector[3] <-  count(df$Sample_time_point == "Batch3")[[2,2]]
  vector[4] <-  count(df$Sample_time_point == "Chemostat1")[[2,2]]
  vector[5] <-  count(df$Sample_time_point == "Chemostat2")[[2,2]]
  vector[6] <-  count(df$Sample_time_point == "Chemostat3")[[2,2]]
  vector[7] <-  count(df$Sample_time_point == "Chemostat4")[[2,2]]
  return(vector)
}


stat_single_compare <- function(df){
  stat_bac<- df %>%
    wilcox_test(Abundance~Sample_time_point,
                p.adjust.method = "fdr",
                paired = FALSE, 
                alternative = "two.sided",
                detailed = TRUE) 
  return(stat_bac)
}




ggplot_boxplot_abundance <- function(df)
{
  ggplot(df, aes(x= Sample_time_point, y= Abundance))+ 
    scale_fill_manual(values = cols_point[2:8]) +
    stat_boxplot(geom ='errorbar', linetype=1, width=0.2) + 
    geom_boxplot(outlier.shape = NA, aes(fill=Sample_time_point)) +
    geom_jitter(show.legend=FALSE, width=0.25, shape=21, fill="black") +
    scale_x_discrete(breaks=group_point[2:8],
                     labels=c(glue("Batch1\nN={counts[1]}"),
                              glue("Batch2\nN={counts[2]}"),
                              glue("Batch3\nN={counts[3]}"),
                              glue("Chemostat1\nN={counts[4]}"),
                              glue("Chemostat2\nN={counts[5]}"),
                              glue("Chemostat3\nN={counts[6]}"),
                              glue("Chemostat4\nN={counts[7]}")
                     )) +
    # scale_y_continuous(limits = c(10, 45))+
    theme_classic() +
    mytheme_alpha_noy
}

sheet <- list()
######################Plot drawing and statistic cal##############################################################################

tab_long_0 <- subset(tab_long, Sample_time_point%in%group_point[2:8])
# 1. HMO####

variable <- "HMO"
df<-subset(tab_long_0, Substrate %in% variable)
stat <- stat_single_compare(df)
stat
counts <- count_number(df)

H_1 <-  ggplot_boxplot_abundance(df)+
  labs(x="",  title=variable)+
  theme(plot.title = element_text(hjust = 0.5))

H_1
sheet[[variable]] <- stat

# 2. Lactose####

variable <- "Lactose"
df<-subset(tab_long_0, Substrate %in% variable)
stat <- stat_single_compare(df)
stat
counts <- count_number(df)

L_1 <-  ggplot_boxplot_abundance(df)+
  labs(x="",  title=variable)+
  theme(plot.title = element_text(hjust = 0.5))

L_1
sheet[[variable]] <- stat





################merge all figures together################

ggarrange(p_inoculum,
          H_1,
          L_1,
          nrow = 1,
          ncol = 3,
          common.legend = T,
          legend = "right",
          widths = c(1.5,5.25,5.25)
) # For export svg 1200x400 size

###################################################################################
# write_xlsx(sheet, "Virome/stat_result/Chemostat/alpha_diversity_substrate_deduct_media.xlsx")

write_xlsx(sheet, "Virome/stat_result/Chemostat/alpha_diversity_substrate.xlsx")
