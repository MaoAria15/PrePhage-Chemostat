setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################Load sequencing data#######################################
##Import files
# source("4_Chemostat.R") # For not deduct media
source("5_Chemostat_deduct_media.R") # For deduct media

ps <- PSV.no.realm
level <- "Family"
# scale <- c(0,100)#origin
scale <- c(0,100) #zoom-in
############################Prepare the taxa#############################

vir.phyl <- tax_glom(ps, level, NArm = FALSE)
ps2 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))
#Create melted dataframe
df <- psmelt(ps2)
#Select last non-empty taxonomic rank
df[df==""] <- NA
df$tax <- apply(df, 1, function(x) tail(na.omit(x), 1))
top<-df %>%
  # filter(level=="Phylum")%>%
  group_by(tax)%>%
  dplyr::summarize(mean_abund=mean(Abundance), .groups = "drop")%>%
  arrange(desc(mean_abund))

#Find top 30 tax
top
top20 <- top$tax[1:6] #Choose number of taxa you want to include in the result

############################Section1: Divided by sample time point#############################
##0.Inoculum##############
Time <- "Inoculum"
ps1<- subset_samples(ps, Sample_time_point %in% Time)

# Get the sample size 
n1 <-  count(ps1@sam_data$Sample_time_point == "Inoculum")[[1,2]]



vir.phyl <- tax_glom(ps1, level, NArm = FALSE)

ps2 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))

ps2 <- merge_samples(ps2, "Sample_time_point")

ps3 <- transform_sample_counts(ps2, function(x) x / sum(x))
#Create melted dataframe
df <- psmelt(ps3)

#Select last non-empty taxonomic rank
df[df==""] <- NA

df$tax <- apply(df, 1, function(x) tail(na.omit(x), 1))

df0 <- df %>%
  mutate(tax = fct_other(df$tax, keep=c(as.matrix(top20))))%>%
  arrange(desc(tax))

#Set order for samples
df0$Sample <- factor(df0$Sample, levels = "Inoculum")


df0$tax <- factor(df0$tax,level=c(top20,"Other"))
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
df0$Abundance <- df0$Abundance*100
p_Inoculum <- ggplot(df0, aes(Sample, Abundance, fill = tax)) + 
  geom_col(width = 0.8) +
  theme_classic() +
  scale_fill_jco(name = "Taxonomy") +
  scale_fill_manual(values=rep(col_vector,10),name=level)+
  scale_x_discrete(labels=c(glue("Inoculum\n(N={n1})"))
  )+
  coord_cartesian(ylim = scale) +
  mytheme_abundance+
  labs(x= "",y="Mean Relative abundance (%)",  title=Time)+
  theme(plot.title = element_text(hjust = 0.5))

p_Inoculum


##1.Batch1##############
Time <- "Batch1"
ps1<- subset_samples(ps, Sample_time_point %in% Time)

# Get the sample size 
n1 <-  count(ps1@sam_data$Substrate == "HMO")[[2,2]]
n2 <-  count(ps1@sam_data$Substrate == "Lactose")[[2,2]]


vir.phyl <- tax_glom(ps1, level, NArm = FALSE)

ps2 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))

ps2 <- merge_samples(ps2, "Substrate")

ps3 <- transform_sample_counts(ps2, function(x) x / sum(x))
#Create melted dataframe
df <- psmelt(ps3)

#Select last non-empty taxonomic rank
df[df==""] <- NA

df$tax <- apply(df, 1, function(x) tail(na.omit(x), 1))

df0 <- df %>%
  mutate(tax = fct_other(df$tax, keep=c(as.matrix(top20))))%>%
  arrange(desc(tax))

#Set order for samples
df0$Sample <- factor(df0$Sample, levels = group_subsrate)


df0$tax <- factor(df0$tax,level=c(top20,"Other"))
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
df0$Abundance <- df0$Abundance*100
B_1 <- ggplot(df0, aes(Sample, Abundance, fill = tax)) + 
  geom_col(width = 0.8) +
  theme_classic() +
  scale_fill_jco(name = "Taxonomy") +
  scale_fill_manual(values=rep(col_vector,10),name=level)+
  scale_x_discrete(labels=c(glue("HMO\n(N={n1})"),
                            glue("Lactose\n(N={n1})"))
  )+
  coord_cartesian(ylim = scale) +
  mytheme_abundance_noy+
  labs(x= "",y="Mean Relative abundance (%)",  title=Time)+
  theme(plot.title = element_text(hjust = 0.5))

B_1

##2.Batch2##############
Time <- "Batch2"
ps1<- subset_samples(ps, Sample_time_point %in% Time)

# Get the sample size 
n1 <-  count(ps1@sam_data$Substrate == "HMO")[[2,2]]
n2 <-  count(ps1@sam_data$Substrate == "Lactose")[[2,2]]


vir.phyl <- tax_glom(ps1, level, NArm = FALSE)

ps2 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))

ps2 <- merge_samples(ps2, "Substrate")

ps3 <- transform_sample_counts(ps2, function(x) x / sum(x))
#Create melted dataframe
df <- psmelt(ps3)

#Select last non-empty taxonomic rank
df[df==""] <- NA

df$tax <- apply(df, 1, function(x) tail(na.omit(x), 1))

df0 <- df %>%
  mutate(tax = fct_other(df$tax, keep=c(as.matrix(top20))))%>%
  arrange(desc(tax))

#Set order for samples
df0$Sample <- factor(df0$Sample, levels = group_subsrate)


df0$tax <- factor(df0$tax,level=c(top20,"Other"))
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
df0$Abundance <- df0$Abundance*100
B_2 <- ggplot(df0, aes(Sample, Abundance, fill = tax)) + 
  geom_col(width = 0.8) +
  theme_classic() +
  scale_fill_jco(name = "Taxonomy") +
  scale_fill_manual(values=rep(col_vector,10),name=level)+
  scale_x_discrete(labels=c(glue("HMO\n(N={n1})"),
                            glue("Lactose\n(N={n1})"))
  )+
  coord_cartesian(ylim = scale) +
  mytheme_abundance_noy+
  labs(x= "",y="Mean Relative abundance (%)",  title=Time)+
  theme(plot.title = element_text(hjust = 0.5))

B_2

##3.Batch3##############
Time <- "Batch3"
ps1<- subset_samples(ps, Sample_time_point %in% Time)

# Get the sample size 
n1 <-  count(ps1@sam_data$Substrate == "HMO")[[2,2]]
n2 <-  count(ps1@sam_data$Substrate == "Lactose")[[2,2]]


vir.phyl <- tax_glom(ps1, level, NArm = FALSE)

ps2 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))

ps2 <- merge_samples(ps2, "Substrate")

ps3 <- transform_sample_counts(ps2, function(x) x / sum(x))
#Create melted dataframe
df <- psmelt(ps3)

#Select last non-empty taxonomic rank
df[df==""] <- NA

df$tax <- apply(df, 1, function(x) tail(na.omit(x), 1))

df0 <- df %>%
  mutate(tax = fct_other(df$tax, keep=c(as.matrix(top20))))%>%
  arrange(desc(tax))

#Set order for samples
df0$Sample <- factor(df0$Sample, levels = group_subsrate)


df0$tax <- factor(df0$tax,level=c(top20,"Other"))
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
df0$Abundance <- df0$Abundance*100
B_3 <- ggplot(df0, aes(Sample, Abundance, fill = tax)) + 
  geom_col(width = 0.8) +
  theme_classic() +
  scale_fill_jco(name = "Taxonomy") +
  scale_fill_manual(values=rep(col_vector,10),name=level)+
  scale_x_discrete(labels=c(glue("HMO\n(N={n1})"),
                            glue("Lactose\n(N={n1})"))
  )+
  coord_cartesian(ylim = scale) +
  mytheme_abundance_noy+
  labs(x= "",y="Mean Relative abundance (%)",  title=Time)+
  theme(plot.title = element_text(hjust = 0.5))

B_3


##4.Chemostat1##############
Time <- "Chemostat1"
ps1<- subset_samples(ps, Sample_time_point %in% Time)

# Get the sample size 
n1 <-  count(ps1@sam_data$Substrate == "HMO")[[2,2]]
n2 <-  count(ps1@sam_data$Substrate == "Lactose")[[2,2]]


vir.phyl <- tax_glom(ps1, level, NArm = FALSE)

ps2 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))

ps2 <- merge_samples(ps2, "Substrate")

ps3 <- transform_sample_counts(ps2, function(x) x / sum(x))
#Create melted dataframe
df <- psmelt(ps3)

#Select last non-empty taxonomic rank
df[df==""] <- NA

df$tax <- apply(df, 1, function(x) tail(na.omit(x), 1))

df0 <- df %>%
  mutate(tax = fct_other(df$tax, keep=c(as.matrix(top20))))%>%
  arrange(desc(tax))

#Set order for samples
df0$Sample <- factor(df0$Sample, levels = group_subsrate)


df0$tax <- factor(df0$tax,level=c(top20,"Other"))
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
df0$Abundance <- df0$Abundance*100
C_1 <- ggplot(df0, aes(Sample, Abundance, fill = tax)) + 
  geom_col(width = 0.8) +
  theme_classic() +
  scale_fill_jco(name = "Taxonomy") +
  scale_fill_manual(values=rep(col_vector,10),name=level)+
  scale_x_discrete(labels=c(glue("HMO\n(N={n1})"),
                            glue("Lactose\n(N={n1})"))
  )+
  coord_cartesian(ylim = scale) +
  mytheme_abundance_noy+
  labs(x= "",y="Mean Relative abundance (%)",  title=Time)+
  theme(plot.title = element_text(hjust = 0.5))

C_1


##5.Chemostat2##############
Time <- "Chemostat2"
ps1<- subset_samples(ps, Sample_time_point %in% Time)

# Get the sample size 
n1 <-  count(ps1@sam_data$Substrate == "HMO")[[2,2]]
n2 <-  count(ps1@sam_data$Substrate == "Lactose")[[2,2]]


vir.phyl <- tax_glom(ps1, level, NArm = FALSE)

ps2 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))

ps2 <- merge_samples(ps2, "Substrate")

ps3 <- transform_sample_counts(ps2, function(x) x / sum(x))
#Create melted dataframe
df <- psmelt(ps3)

#Select last non-empty taxonomic rank
df[df==""] <- NA

df$tax <- apply(df, 1, function(x) tail(na.omit(x), 1))

df0 <- df %>%
  mutate(tax = fct_other(df$tax, keep=c(as.matrix(top20))))%>%
  arrange(desc(tax))

#Set order for samples
df0$Sample <- factor(df0$Sample, levels = group_subsrate)


df0$tax <- factor(df0$tax,level=c(top20,"Other"))
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
df0$Abundance <- df0$Abundance*100
C_2 <- ggplot(df0, aes(Sample, Abundance, fill = tax)) + 
  geom_col(width = 0.8) +
  theme_classic() +
  scale_fill_jco(name = "Taxonomy") +
  scale_fill_manual(values=rep(col_vector,10),name=level)+
  scale_x_discrete(labels=c(glue("HMO\n(N={n1})"),
                            glue("Lactose\n(N={n1})"))
  )+
  coord_cartesian(ylim = scale) +
  mytheme_abundance_noy+
  labs(x= "",y="Mean Relative abundance (%)",  title=Time)+
  theme(plot.title = element_text(hjust = 0.5))

C_2

##6.Chemostat3##############
Time <- "Chemostat3"
ps1<- subset_samples(ps, Sample_time_point %in% Time)

# Get the sample size 
n1 <-  count(ps1@sam_data$Substrate == "HMO")[[2,2]]
n2 <-  count(ps1@sam_data$Substrate == "Lactose")[[2,2]]


vir.phyl <- tax_glom(ps1, level, NArm = FALSE)

ps2 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))

ps2 <- merge_samples(ps2, "Substrate")

ps3 <- transform_sample_counts(ps2, function(x) x / sum(x))
#Create melted dataframe
df <- psmelt(ps3)

#Select last non-empty taxonomic rank
df[df==""] <- NA

df$tax <- apply(df, 1, function(x) tail(na.omit(x), 1))

df0 <- df %>%
  mutate(tax = fct_other(df$tax, keep=c(as.matrix(top20))))%>%
  arrange(desc(tax))

#Set order for samples
df0$Sample <- factor(df0$Sample, levels = group_subsrate)


df0$tax <- factor(df0$tax,level=c(top20,"Other"))
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
df0$Abundance <- df0$Abundance*100
C_3 <- ggplot(df0, aes(Sample, Abundance, fill = tax)) + 
  geom_col(width = 0.8) +
  theme_classic() +
  scale_fill_jco(name = "Taxonomy") +
  scale_fill_manual(values=rep(col_vector,10),name=level)+
  scale_x_discrete(labels=c(glue("HMO\n(N={n1})"),
                            glue("Lactose\n(N={n1})"))
  )+
  coord_cartesian(ylim = scale) +
  mytheme_abundance_noy+
  labs(x= "",y="Mean Relative abundance (%)",  title=Time)+
  theme(plot.title = element_text(hjust = 0.5))

C_3


##7.Chemostat4##############
Time <- "Chemostat4"
ps1<- subset_samples(ps, Sample_time_point %in% Time)

# Get the sample size 
n1 <-  count(ps1@sam_data$Substrate == "HMO")[[2,2]]
n2 <-  count(ps1@sam_data$Substrate == "Lactose")[[2,2]]


vir.phyl <- tax_glom(ps1, level, NArm = FALSE)

ps2 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))

ps2 <- merge_samples(ps2, "Substrate")

ps3 <- transform_sample_counts(ps2, function(x) x / sum(x))
#Create melted dataframe
df <- psmelt(ps3)

#Select last non-empty taxonomic rank
df[df==""] <- NA

df$tax <- apply(df, 1, function(x) tail(na.omit(x), 1))

df0 <- df %>%
  mutate(tax = fct_other(df$tax, keep=c(as.matrix(top20))))%>%
  arrange(desc(tax))

#Set order for samples
df0$Sample <- factor(df0$Sample, levels = group_subsrate)


df0$tax <- factor(df0$tax,level=c(top20,"Other"))
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
df0$Abundance <- df0$Abundance*100
C_4 <- ggplot(df0, aes(Sample, Abundance, fill = tax)) + 
  geom_col(width = 0.8) +
  theme_classic() +
  scale_fill_jco(name = "Taxonomy") +
  scale_fill_manual(values=rep(col_vector,10),name=level)+
  scale_x_discrete(labels=c(glue("HMO\n(N={n1})"),
                            glue("Lactose\n(N={n1})"))
  )+
  coord_cartesian(ylim = scale) +
  mytheme_abundance_noy+
  labs(x= "",y="Mean Relative abundance (%)",  title=Time)+
  theme(plot.title = element_text(hjust = 0.5))

C_4
##merge#####
ggarrange(p_Inoculum,
          B_1,B_2,B_3,
          C_1,C_2,C_3,C_4,
          nrow = 1,
          ncol = 8,
          common.legend = T,
          legend = "right",
          widths = c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5)) #Size:1000*400



############################Section2: Divided by Substrate#############################

##1.HMO##############
Sub<- "HMO"
ps1<- subset_samples(ps, Substrate %in% Sub)
ps1<- subset_samples(ps1, Sample_time_point %in% group_point[2:8])
# Get the sample size 
n1 <-  count(ps1@sam_data$Sample_time_point == "Batch1")[[2,2]]
n2 <-  count(ps1@sam_data$Sample_time_point == "Batch2")[[2,2]]
n3 <-  count(ps1@sam_data$Sample_time_point == "Batch3")[[2,2]]
n4 <-  count(ps1@sam_data$Sample_time_point == "Chemostat1")[[2,2]]
n5 <-  count(ps1@sam_data$Sample_time_point == "Chemostat2")[[2,2]]
n6 <-  count(ps1@sam_data$Sample_time_point == "Chemostat3")[[2,2]]
n7 <-  count(ps1@sam_data$Sample_time_point == "Chemostat4")[[2,2]]


vir.phyl <- tax_glom(ps1, level, NArm = FALSE)

ps2 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))

ps2 <- merge_samples(ps2, "Sample_time_point")

ps3 <- transform_sample_counts(ps2, function(x) x / sum(x))
#Create melted dataframe
df <- psmelt(ps3)

#Select last non-empty taxonomic rank
df[df==""] <- NA

df$tax <- apply(df, 1, function(x) tail(na.omit(x), 1))

df0 <- df %>%
  mutate(tax = fct_other(df$tax, keep=c(as.matrix(top20))))%>%
  arrange(desc(tax))

#Set order for samples
df0$Sample <- factor(df0$Sample, levels = group_point)


df0$tax <- factor(df0$tax,level=c(top20,"Other"))
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
df0$Abundance <- df0$Abundance*100
p_HMO <- ggplot(df0, aes(Sample, Abundance, fill = tax)) + 
  geom_col(width = 0.8) +
  theme_classic() +
  scale_fill_jco(name = "Taxonomy") +
  scale_fill_manual(values=rep(col_vector,10),name=level)+
  scale_x_discrete(labels=c(glue("Batch1\nN={n1}"),
                            glue("Batch2\nN={n2}"),
                            glue("Batch3\nN={n3}"),
                            glue("Chemostat1\nN={n4}"),
                            glue("Chemostat2\nN={n5}"),
                            glue("Chemostat3\nN={n6}"),
                            glue("Chemostat4\nN={n7}"))
  )+
  coord_cartesian(ylim = scale) +
  mytheme_abundance_noy+
  labs(x= "",y="Mean Relative abundance (%)",  title=Sub)+
  theme(plot.title = element_text(hjust = 0.5))

p_HMO
##2.Lactose#############
Sub<- "Lactose"
ps1<- subset_samples(ps, Substrate %in% Sub)
ps1<- subset_samples(ps1, Sample_time_point %in% group_point[2:8])

# Get the sample size 
n1 <-  count(ps1@sam_data$Sample_time_point == "Batch1")[[2,2]]
n2 <-  count(ps1@sam_data$Sample_time_point == "Batch2")[[2,2]]
n3 <-  count(ps1@sam_data$Sample_time_point == "Batch3")[[2,2]]
n4 <-  count(ps1@sam_data$Sample_time_point == "Chemostat1")[[2,2]]
n5 <-  count(ps1@sam_data$Sample_time_point == "Chemostat2")[[2,2]]
n6 <-  count(ps1@sam_data$Sample_time_point == "Chemostat3")[[2,2]]
n7 <-  count(ps1@sam_data$Sample_time_point == "Chemostat4")[[2,2]]


vir.phyl <- tax_glom(ps1, level, NArm = FALSE)

ps2 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))

ps2 <- merge_samples(ps2, "Sample_time_point")

ps3 <- transform_sample_counts(ps2, function(x) x / sum(x))
#Create melted dataframe
df <- psmelt(ps3)

#Select last non-empty taxonomic rank
df[df==""] <- NA

df$tax <- apply(df, 1, function(x) tail(na.omit(x), 1))

df0 <- df %>%
  mutate(tax = fct_other(df$tax, keep=c(as.matrix(top20))))%>%
  arrange(desc(tax))

#Set order for samples
df0$Sample <- factor(df0$Sample, levels = group_point)


df0$tax <- factor(df0$tax,level=c(top20,"Other"))
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
df0$Abundance <- df0$Abundance*100
p_Lactose <- ggplot(df0, aes(Sample, Abundance, fill = tax)) + 
  geom_col(width = 0.8) +
  theme_classic() +
  scale_fill_jco(name = "Taxonomy") +
  scale_fill_manual(values=rep(col_vector,10),name=level)+
  scale_x_discrete(labels=c(glue("Batch1\nN={n1}"),
                            glue("Batch2\nN={n2}"),
                            glue("Batch3\nN={n3}"),
                            glue("Chemostat1\nN={n4}"),
                            glue("Chemostat2\nN={n5}"),
                            glue("Chemostat3\nN={n6}"),
                            glue("Chemostat4\nN={n7}"))
  )+
  coord_cartesian(ylim = scale) +
  mytheme_abundance_noy+
  labs(x= "",y="Mean Relative abundance (%)",  title=Sub)+
  theme(plot.title = element_text(hjust = 0.5))

p_Lactose

##merge#####
ggarrange(p_Inoculum,
          p_HMO,
          p_Lactose,
          nrow = 1,
          ncol = 3,
          common.legend = T,
          legend = "right",
          widths = c(1.5,5.25,5.25)) #Size:1200*400

############################Section3: Divided by Substrate by samples#############################

##1.HMO##############
Sub<- "HMO"

ps1<- subset_samples(ps, Substrate %in% Sub)

vir.phyl <- tax_glom(ps1, level, NArm = FALSE)

ps2 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))

#Create melted dataframe
df <- psmelt(ps2)

#Select last non-empty taxonomic rank
df[df==""] <- NA

df$tax <- apply(df, 1, function(x) tail(na.omit(x), 1))

df0 <- df %>%
  mutate(tax = fct_other(df$tax, keep=c(as.matrix(top20))))%>%
  arrange(desc(tax))
df0$tax <- factor(df0$tax,level=c(top20,"Other"))

#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
p_HMO <- ggplot(df0, aes(Propagation , Abundance, fill = tax)) + 
  geom_col(width = 0.8) +
  theme_classic() +
  scale_fill_jco(name = "Taxonomy") +
  scale_fill_manual(values=rep(col_vector,10),name=level)+
  coord_cartesian(ylim = scale) +
  mytheme_abundance_yvertival+
  labs(x= "",y="Mean Relative abundance (%)",  title=Sub)+
  theme(plot.title = element_text(hjust = 0.5))

p_HMO

##2.Lactose#############
Sub<- "Lactose"
ps1<- subset_samples(ps, Substrate %in% Sub)

vir.phyl <- tax_glom(ps1, level, NArm = FALSE)

ps2 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))

#Create melted dataframe
df <- psmelt(ps2)

#Select last non-empty taxonomic rank
df[df==""] <- NA

df$tax <- apply(df, 1, function(x) tail(na.omit(x), 1))

df0 <- df %>%
  mutate(tax = fct_other(df$tax, keep=c(as.matrix(top20))))%>%
  arrange(desc(tax))

df0$tax <- factor(df0$tax,level=c(top20,"Other"))
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
p_Lactose <- ggplot(df0, aes(Propagation , Abundance, fill = tax)) + 
  geom_col(width = 0.8) +
  theme_classic() +
  scale_fill_jco(name = "Taxonomy") +
  scale_fill_manual(values=rep(col_vector,10),name=level)+
  coord_cartesian(ylim = scale) +
  mytheme_abundance_yvertival+
  labs(x= "",y="Mean Relative abundance (%)",  title=Sub)+
  theme(plot.title = element_text(hjust = 0.5))

p_Lactose

##merge#####
ggarrange(p_HMO,
          p_Lactose,
          nrow = 1,
          ncol = 2,
          common.legend = T,
          legend = "right",
          widths = c(10,5)) #Size:1500*400

