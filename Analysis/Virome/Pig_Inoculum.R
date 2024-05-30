setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################Load sequencing data#######################################

##Import files



source("4_Donors(chemostat)_Inoculum(pig).R")


ps <- subset_samples(PSV.no.realm, Sample_origin %in% "Pig_inoculum")

###########################1. alpha diversity#####################################
alpha_met <- c("Observed") 


standf = function(x, t=total) round(t * (x / sum(x)))
total = median(sample_sums(ps))
PSV.R = transform_sample_counts(ps, standf)

rich <- estimate_richness(PSV.R, measures = alpha_met)
tab <- subset(rich, select = alpha_met) 
index <- match(rownames(rich), rownames(sample_data(PSV.R)))
tab$Group <-sample_data(PSV.R)$Group[index]


tab$Group <- factor(tab$Group,levels = c("Inoculum_FVT","Inoculum_CVT","Inoculum_CVT_MO","Inoculum_CON"))


n1 <-  count(tab$Group == "Inoculum_FVT")[[2,2]]
n2 <-  count(tab$Group == "Inoculum_CVT")[[2,2]]
n3 <-  count(tab$Group == "Inoculum_CVT_MO")[[2,2]]
n4 <-  count(tab$Group == "Inoculum_CON")[[2,2]]

p_alpha <- ggplot(tab, aes(x = Group, y = Observed))+
  geom_bar(stat = "identity", position = "dodge",aes(fill=Group),color="black") +
  scale_fill_manual(values = simone_cols) +
  labs(x="", y="Observed alpha diversity", title="Virome") +
  scale_x_discrete(labels=c(glue("Inoculum\nFVT\nN={n1}"),
                            glue("Inoculum\nCVT\nN={n2}"),
                            glue("Inoculum\nCVT_MO\nN={n3}"),
                            glue("Inoculum\nCON\nN={n4}"))
  )+
  scale_y_continuous(limits = c(0,450))+
  theme_classic() +
  geom_text(aes(label=Observed), vjust=-0.5,size=3)+
  mytheme_beta

p_alpha
#############################2. abundant barplot##################################

level <- "Family"

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
top20 <- top$tax[1:5] #Choose number of taxa you want to include in the result


df0 <- df %>%
  mutate(tax = fct_other(df$tax, keep=c(as.matrix(top20))))%>%
  arrange(desc(tax))

#Set order for samples
df0$Sample <- factor(df0$Group, levels = c("Inoculum_FVT","Inoculum_CVT","Inoculum_CVT_MO","Inoculum_CON"))


df0$tax <- factor(df0$tax,level=c(top20,"Other"))
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
p_barplot <- ggplot(df0, aes(Sample, Abundance, fill = tax)) + 
  geom_col(width = 0.8) +
  theme_classic() +
  scale_fill_jco(name = "Taxonomy") +
  scale_fill_manual(values=rep(col_vector,10),name=level)+
  scale_x_discrete(labels=c(glue("Inoculum\nFVT\nN={n1}"),
                            glue("Inoculum\nCVT\nN={n2}"),
                            glue("Inoculum\nCVT_MO\nN={n3}"),
                            glue("Inoculum\nCON\nN={n4}"))
  )+
  coord_cartesian(ylim = c(0, 100)) +
  mytheme_abundance+
  labs(x= "",y="Mean Relative abundance (%)",  title="Virome abundance")+
  theme(plot.title = element_text(hjust = 0.5))

p_barplot


#############################3. Host abundant barplot##################################
ps <- subset_samples(PSV.HOST, Sample_origin %in% "Pig_inoculum")
level <- "Genus"

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
top20 <- top$tax[1:15] #Choose number of taxa you want to include in the result


df0 <- df %>%
  mutate(tax = fct_other(df$tax, keep=c(as.matrix(top20))))%>%
  arrange(desc(tax))

#Set order for samples
df0$Sample <- factor(df0$Group, levels = c("Inoculum_FVT","Inoculum_CVT","Inoculum_CVT_MO","Inoculum_CON"))


df0$tax <- factor(df0$tax,level=c(top20,"Other"))
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
p_barplot <- ggplot(df0, aes(Sample, Abundance, fill = tax)) + 
  geom_col(width = 0.8) +
  theme_classic() +
  scale_fill_jco(name = "Taxonomy") +
  scale_fill_manual(values=rep(col_vector,10),name=level)+
  scale_x_discrete(labels=c(glue("Inoculum\nFVT\nN={n1}"),
                            glue("Inoculum\nCVT\nN={n2}"),
                            glue("Inoculum\nCVT_MO\nN={n3}"),
                            glue("Inoculum\nCON\nN={n4}"))
  )+
  coord_cartesian(ylim = c(0, 100)) +
  mytheme_abundance+
  labs(x= "",y="Mean Relative abundance (%)",  title="Host prediction")+
  theme(plot.title = element_text(hjust = 0.5))

p_barplot
