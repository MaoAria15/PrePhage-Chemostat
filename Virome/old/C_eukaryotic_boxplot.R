
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("4_Chemostat.R")


ps1 <- subset_samples(PSV.no.realm, Sample_time_point %in% group_point) 


ps2<- ps1

  taxonomy <- read.csv('Data/Virome/eukaryotic_taxonomy.csv' ,row.names = 1)
tax_table(ps2) <- as.matrix(taxonomy)


sum_counts1  <- colSums(otu_table(ps1))
sum_counts2  <- colSums(otu_table(ps2))
percentage <- as.data.frame((sum_counts2 / sum_counts1) * 100)

mapping <- as.data.frame(sample_data(ps1))
percentage$Eukaryotic_Percentage<- percentage$`(sum_counts2/sum_counts1) * 100`
percentage$Sample_time_point<- mapping$Sample_time_point
percentage$Substrate<- mapping$Substrate

sheet <- list()
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
    wilcox_test(Eukaryotic_Percentage~Substrate,
                # p.adjust.method = "fdr",
                paired = FALSE, 
                alternative = "two.sided",
                detailed = TRUE) 
  return(stat_bac)
}



ggplot_boxplot_abundance <- function(df)
{
  ggplot(df, aes(x= Substrate, y= Eukaryotic_Percentage))+ 
    scale_fill_manual(values = col_substrate) +
    stat_boxplot(geom ='errorbar', linetype=1, width=0.2) + 
    geom_boxplot(outlier.shape = NA, aes(fill=Substrate)) +
    geom_jitter(show.legend=FALSE, width=0.25, shape=21, fill="black") +
    scale_x_discrete(breaks=group_subsrate,
                     labels=c(glue("HMO\nN={counts[1]}"),
                              glue("Lactose\nN={counts[2]}"))) +
    scale_y_continuous(limits = c(0, 1))+
    theme_classic() +
    mytheme_alpha_noy
}

sheet <- list()
######################Plot drawing and statistic cal##############################################################################
# 0. Inoculum####

variable <- "Inoculum"
df<-subset(percentage, Sample_time_point %in% variable)
counts <- count(df$Sample_time_point == variable)[[1,2]]

p_inoculum  <-  ggplot(df, aes(x= Sample_time_point, y= Eukaryotic_Percentage))+ 
  stat_boxplot(geom ='errorbar', linetype=1, width=0.2) + 
  geom_boxplot(outlier.shape = NA, aes(fill=Substrate)) +
  scale_fill_manual(values = cols_point[1]) + 
  geom_jitter(show.legend=FALSE, width=0.25, shape=21, fill="black") +
  scale_x_discrete(labels=c(glue("Inoculum\nN={counts[1]}"))) +
  theme_classic() +
  mytheme_beta+
  labs(x="", y="Eukaryotic viruses abundance (%)", title=variable)+
  scale_y_continuous(limits = c(0, 1))+
  theme(plot.title = element_text(hjust = 0.5))

p_inoculum

# 1. Batch1####

variable <- "Batch1"
df<-subset(percentage, Sample_time_point %in% variable)
stat <- stat_single_compare(df)
stat
counts <- count_number(df)

B_1 <-  ggplot_boxplot_abundance(df)+
  labs(x="",  title=variable)+
  theme(plot.title = element_text(hjust = 0.5))

B_1


# 2. Batch2####

variable <- "Batch2"
df<-subset(percentage, Sample_time_point %in% variable)
stat <- stat_single_compare(df)
stat
counts <- count_number(df)

B_2 <-  ggplot_boxplot_abundance(df)+
  labs(x="",  title=variable)+
  theme(plot.title = element_text(hjust = 0.5))

B_2



# 3. Batch3####

variable <- "Batch3"
df<-subset(percentage, Sample_time_point %in% variable)
stat <- stat_single_compare(df)
stat
counts <- count_number(df)

B_3 <-  ggplot_boxplot_abundance(df)+
  labs(x="",  title=variable)+
  theme(plot.title = element_text(hjust = 0.5))

B_3



# 4. Chemostat1####

variable <- "Chemostat1"
df<-subset(percentage, Sample_time_point %in% variable)
stat <- stat_single_compare(df)
stat
counts <- count_number(df)

C_1 <-  ggplot_boxplot_abundance(df)+
  labs(x="",  title=variable)+
  theme(plot.title = element_text(hjust = 0.5))

C_1



# 5. Chemostat2####

variable <- "Chemostat2"
df<-subset(percentage, Sample_time_point %in% variable)
stat <- stat_single_compare(df)
stat
counts <- count_number(df)

C_2 <-  ggplot_boxplot_abundance(df)+
  labs(x="",  title=variable)+
  theme(plot.title = element_text(hjust = 0.5))

C_2


# 6. Chemostat3####

variable <- "Chemostat3"
df<-subset(percentage, Sample_time_point %in% variable)
stat <- stat_single_compare(df)
stat
counts <- count_number(df)

C_3 <-  ggplot_boxplot_abundance(df)+
  labs(x="",  title=variable)+
  theme(plot.title = element_text(hjust = 0.5))

C_3


# 7. Chemostat4####

variable <- "Chemostat4"
df<-subset(percentage, Sample_time_point %in% variable)
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
    wilcox_test(Eukaryotic_Percentage~Sample_time_point,
                p.adjust.method = "fdr",
                paired = FALSE, 
                alternative = "two.sided",
                detailed = TRUE) 
  return(stat_bac)
}




ggplot_boxplot_abundance <- function(df)
{
  ggplot(df, aes(x= Sample_time_point, y= Eukaryotic_Percentage))+ 
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
    scale_y_continuous(limits = c(0, 1))+
    theme_classic() +
    mytheme_alpha_noy
}

sheet <- list()
######################Plot drawing and statistic cal##############################################################################

tab_long_0 <- subset(percentage, Sample_time_point%in%group_point[2:8])
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
write_xlsx(sheet, "Virome/stat_result/Chemostat/eukaryotic_abundance_substrate.xlsx")
