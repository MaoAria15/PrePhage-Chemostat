setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################Load sequencing data#######################################
##Import files
source("4_simone.R")

ps <- subset_samples(PSV.no.realm, Group %in% S_groups)
level <- "Family"
############################Prepare the taxa#############################




# Get the sample size 
n1 <-  count(ps@sam_data$Group == "FVT")[[2,2]]
n2 <-  count(ps@sam_data$Group == "CVT")[[2,2]]
n3 <-  count(ps@sam_data$Group == "CVT_MO")[[2,2]]
n4 <-  count(ps@sam_data$Group == "CON")[[2,2]]


vir.phyl <- tax_glom(ps, level, NArm = FALSE)

ps2 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))

ps2 <- merge_samples(ps2, "Group")

ps3 <- transform_sample_counts(ps2, function(x) x / sum(x))
#Create melted dataframe
df <- psmelt(ps3)

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
df0$Sample <- factor(df0$Sample, levels = S_groups)


df0$tax <- factor(df0$tax,level=c(top20,"Other"))
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
df0$Abundance <- df0$Abundance*100
p_barplot <- ggplot(df0, aes(Sample, Abundance, fill = tax)) + 
  geom_col(width = 0.8) +
  theme_classic() +
  scale_fill_jco(name = "Taxonomy") +
  scale_fill_manual(values=rep(col_vector,10),name=level)+
  scale_x_discrete(labels=c(glue("FVT\nN={n1}"),
                            glue("CVT\nN={n2}"),
                            glue("CVT_MO\nN={n3}"),
                            glue("CON\nN={n4}"))
  )+
  coord_cartesian(ylim = c(0, 100)) +
  mytheme_abundance+
  labs(x= "",y="Mean Relative abundance (%)",  title="Virome")+
  theme(plot.title = element_text(hjust = 0.5))

p_barplot

