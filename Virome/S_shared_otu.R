setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("4_Simone.R")
View(PSV.no.realm@sam_data)
##################################1. FVT############################################

ps_con_inoculum <- subset_samples(PSV.no.realm, Group %in% "Inoculum_FVT")
ps_con <- subset_samples(PSV.no.realm, Group %in% "FVT")



#Subset the Inoculum_CON
otu_sums <- taxa_sums(ps_con_inoculum)
otus_to_keep <- otu_sums >0
ps_con_inoculum <- prune_taxa(otus_to_keep, ps_con_inoculum)


#Subset the CON
otu_sums <- taxa_sums(ps_con)
otus_to_keep <- otu_sums >0
ps_con <- prune_taxa(otus_to_keep, ps_con)

#Get otu names
otu_con_inoculum <- rownames(otu_table(ps_con_inoculum))                                       #
otu_con <- rownames(otu_table(ps_con))

#Find unique otu of CON
share_con <- intersect(otu_con_inoculum,otu_con) 
otu_con <- otu_con[!otu_con %in% share_con]    




# prepare the phyloseq
# Merge the phyloseq objects
ps_merged <- merge_phyloseq(ps_con_inoculum, ps_con)


#get the tax_table
taxonomy <- as.data.frame(tax_table(ps_merged))


# Add a new column 'Source' initialized with 'nono'
taxonomy$Source <- 'None'

# Update 'Source' based on the presence in otu vectors
taxonomy$Source[rownames(taxonomy) %in% otu_con_inoculum] <- 'Inoculum'
taxonomy$Source[rownames(taxonomy) %in% otu_con] <- 'From_pig'


# Convert the modified data frame back to a Taxonomy Table
taxonomy <- tax_table(as.matrix(taxonomy))
# Update the original phyloseq object with the new taxonomy table

tax_table(ps_merged) <- taxonomy

vir.phyl <- tax_glom(ps_merged, "Source", NArm = FALSE)
ps0 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))

ps1 <- merge_samples(ps0, "Group")

ps1 <- transform_sample_counts(ps1, function(x) x / sum(x))


#Create melted dataframe
df <- psmelt(ps1)



#Select last non-empty taxonomic rank
df[df==""] <- NA


#Set order for samples
df$Sample <- factor(df$Sample, levels = c("Inoculum_FVT","FVT"))


df$Abundance <- df$Abundance*100

df$Source <- factor(df$Source,level = c("Inoculum","From_pig"))
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

p_fvt <- ggplot(df, aes(Sample, Abundance, fill = Source)) + 
  geom_col(width = 0.8) +
  theme_classic() +
  scale_fill_jco(name = "Taxonomy") +
  # facet_wrap("Sample_time_point")+
  scale_fill_manual(values=rep(col_vector,10),name="Source")+
  scale_x_discrete(breaks=c("Inoculum_FVT","FVT"),
                   labels=c(glue("Inoculum_FVT"),
                            glue("FVT")))+
  coord_cartesian(ylim = c(0, 100)) +
  mytheme_abundance+
  # theme_classic() + 
  ylab("Mean Relative abundance (%)") +
  xlab("FVT")


p_fvt
##################################2. CVT############################################

ps_con_inoculum <- subset_samples(PSV.no.realm, Group %in% "Inoculum_CVT")
ps_con <- subset_samples(PSV.no.realm, Group %in% "CVT")



#Subset the Inoculum_CON
otu_sums <- taxa_sums(ps_con_inoculum)
otus_to_keep <- otu_sums >0
ps_con_inoculum <- prune_taxa(otus_to_keep, ps_con_inoculum)


#Subset the CON
otu_sums <- taxa_sums(ps_con)
otus_to_keep <- otu_sums >0
ps_con <- prune_taxa(otus_to_keep, ps_con)

#Get otu names
otu_con_inoculum <- rownames(otu_table(ps_con_inoculum))                                       #
otu_con <- rownames(otu_table(ps_con))

#Find unique otu of CON
share_con <- intersect(otu_con_inoculum,otu_con) 
otu_con <- otu_con[!otu_con %in% share_con]    




# prepare the phyloseq
# Merge the phyloseq objects
ps_merged <- merge_phyloseq(ps_con_inoculum, ps_con)


#get the tax_table
taxonomy <- as.data.frame(tax_table(ps_merged))


# Add a new column 'Source' initialized with 'nono'
taxonomy$Source <- 'None'

# Update 'Source' based on the presence in otu vectors
taxonomy$Source[rownames(taxonomy) %in% otu_con_inoculum] <- 'Inoculum'
taxonomy$Source[rownames(taxonomy) %in% otu_con] <- 'From_pig'


# Convert the modified data frame back to a Taxonomy Table
taxonomy <- tax_table(as.matrix(taxonomy))
# Update the original phyloseq object with the new taxonomy table

tax_table(ps_merged) <- taxonomy

vir.phyl <- tax_glom(ps_merged, "Source", NArm = FALSE)
ps0 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))

ps1 <- merge_samples(ps0, "Group")

ps1 <- transform_sample_counts(ps1, function(x) x / sum(x))


#Create melted dataframe
df <- psmelt(ps1)



#Select last non-empty taxonomic rank
df[df==""] <- NA


#Set order for samples
df$Sample <- factor(df$Sample, levels = c("Inoculum_CVT","CVT"))


df$Abundance <- df$Abundance*100

df$Source <- factor(df$Source,level = c("Inoculum","From_pig"))
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

p_cvt <- ggplot(df, aes(Sample, Abundance, fill = Source)) + 
  geom_col(width = 0.8) +
  theme_classic() +
  scale_fill_jco(name = "Taxonomy") +
  # facet_wrap("Sample_time_point")+
  scale_fill_manual(values=rep(col_vector,10),name="Source")+
  scale_x_discrete(breaks=c("Inoculum_CVT","CVT"),
                   labels=c(glue("Inoculum_CVT"),
                            glue("CVT")))+
  coord_cartesian(ylim = c(0, 100)) +
  mytheme_abundance+
  # theme_classic() + 
  ylab("Mean Relative abundance (%)") +
  xlab("CVT")


p_cvt

##################################3. CVT-MO############################################

ps_con_inoculum <- subset_samples(PSV.no.realm, Group %in% "Inoculum_CVT_MO")
ps_con <- subset_samples(PSV.no.realm, Group %in% "CVT_MO")



#Subset the Inoculum_CON
otu_sums <- taxa_sums(ps_con_inoculum)
otus_to_keep <- otu_sums >0
ps_con_inoculum <- prune_taxa(otus_to_keep, ps_con_inoculum)


#Subset the CON
otu_sums <- taxa_sums(ps_con)
otus_to_keep <- otu_sums >0
ps_con <- prune_taxa(otus_to_keep, ps_con)

#Get otu names
otu_con_inoculum <- rownames(otu_table(ps_con_inoculum))                                       #
otu_con <- rownames(otu_table(ps_con))

#Find unique otu of CON
share_con <- intersect(otu_con_inoculum,otu_con) 
otu_con <- otu_con[!otu_con %in% share_con]    




# prepare the phyloseq
# Merge the phyloseq objects
ps_merged <- merge_phyloseq(ps_con_inoculum, ps_con)


#get the tax_table
taxonomy <- as.data.frame(tax_table(ps_merged))


# Add a new column 'Source' initialized with 'nono'
taxonomy$Source <- 'None'

# Update 'Source' based on the presence in otu vectors
taxonomy$Source[rownames(taxonomy) %in% otu_con_inoculum] <- 'Inoculum'
taxonomy$Source[rownames(taxonomy) %in% otu_con] <- 'From_pig'


# Convert the modified data frame back to a Taxonomy Table
taxonomy <- tax_table(as.matrix(taxonomy))
# Update the original phyloseq object with the new taxonomy table

tax_table(ps_merged) <- taxonomy

vir.phyl <- tax_glom(ps_merged, "Source", NArm = FALSE)
ps0 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))

ps1 <- merge_samples(ps0, "Group")

ps1 <- transform_sample_counts(ps1, function(x) x / sum(x))


#Create melted dataframe
df <- psmelt(ps1)



#Select last non-empty taxonomic rank
df[df==""] <- NA


#Set order for samples
df$Sample <- factor(df$Sample, levels = c("Inoculum_CVT_MO","CVT_MO"))


df$Abundance <- df$Abundance*100

df$Source <- factor(df$Source,level = c("Inoculum","From_pig"))
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

p_CVTMO <- ggplot(df, aes(Sample, Abundance, fill = Source)) + 
  geom_col(width = 0.8) +
  theme_classic() +
  scale_fill_jco(name = "Taxonomy") +
  # facet_wrap("Sample_time_point")+
  scale_fill_manual(values=rep(col_vector,10),name="Source")+
  scale_x_discrete(breaks=c("Inoculum_CVT_MO","CVT_MO"),
                   labels=c(glue("Inoculum_CVT_MO"),
                            glue("CVT_MO")))+
  coord_cartesian(ylim = c(0, 100)) +
  mytheme_abundance+
  # theme_classic() + 
  ylab("Mean Relative abundance (%)") +
  xlab("CVT_MO")


p_CVTMO
##################################4. CON############################################

ps_con_inoculum <- subset_samples(PSV.no.realm, Group %in% "Inoculum_CON")
ps_con <- subset_samples(PSV.no.realm, Group %in% "CON")



#Subset the Inoculum_CON
otu_sums <- taxa_sums(ps_con_inoculum)
otus_to_keep <- otu_sums >0
ps_con_inoculum <- prune_taxa(otus_to_keep, ps_con_inoculum)


#Subset the CON
otu_sums <- taxa_sums(ps_con)
otus_to_keep <- otu_sums >0
ps_con <- prune_taxa(otus_to_keep, ps_con)

#Get otu names
otu_con_inoculum <- rownames(otu_table(ps_con_inoculum))                                       #
otu_con <- rownames(otu_table(ps_con))

#Find unique otu of CON
share_con <- intersect(otu_con_inoculum,otu_con) 
otu_con <- otu_con[!otu_con %in% share_con]    




# prepare the phyloseq
# Merge the phyloseq objects
ps_merged <- merge_phyloseq(ps_con_inoculum, ps_con)


#get the tax_table
taxonomy <- as.data.frame(tax_table(ps_merged))


# Add a new column 'Source' initialized with 'nono'
taxonomy$Source <- 'None'

# Update 'Source' based on the presence in otu vectors
taxonomy$Source[rownames(taxonomy) %in% otu_con_inoculum] <- 'Inoculum'
taxonomy$Source[rownames(taxonomy) %in% otu_con] <- 'From_pig'


# Convert the modified data frame back to a Taxonomy Table
taxonomy <- tax_table(as.matrix(taxonomy))
# Update the original phyloseq object with the new taxonomy table

tax_table(ps_merged) <- taxonomy

vir.phyl <- tax_glom(ps_merged, "Source", NArm = FALSE)
ps0 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))

ps1 <- merge_samples(ps0, "Group")

ps1 <- transform_sample_counts(ps1, function(x) x / sum(x))


#Create melted dataframe
df <- psmelt(ps1)



#Select last non-empty taxonomic rank
df[df==""] <- NA


#Set order for samples
df$Sample <- factor(df$Sample, levels = c("Inoculum_CON","CON"))


df$Abundance <- df$Abundance*100

df$Source <- factor(df$Source,level = c("Inoculum","From_pig"))
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

p_con <- ggplot(df, aes(Sample, Abundance, fill = Source)) + 
  geom_col(width = 0.8) +
  theme_classic() +
  scale_fill_jco(name = "Taxonomy") +
  # facet_wrap("Sample_time_point")+
  scale_fill_manual(values=rep(col_vector,10),name="Source")+
  scale_x_discrete(breaks=c("Inoculum_CON","CON"),
                   labels=c(glue("Inoculum_CON"),
                            glue("CON")))+
  coord_cartesian(ylim = c(0, 100)) +
  mytheme_abundance+
  # theme_classic() + 
  ylab("Mean Relative abundance (%)") +
  xlab("CON")


p_con



ggarrange(p_fvt,
          p_cvt,
          p_CVTMO,
          p_con,
          nrow = 1,
          ncol = 4,
          common.legend = T)

