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

tax <- data.frame(Group=rep("Viruse",length(taxonomy$Source)),
                  Source=taxonomy$Source)

rownames(tax) <- rownames(taxonomy)




# Convert the modified data frame back to a Taxonomy Table
tax <- tax_table(as.matrix(tax))
# Update the original phyloseq object with the new taxonomy table

tax_table(ps_merged) <- tax

ps_merged <- subset_samples(ps_merged, Sample_origin %in% "Pig")

vir.phyl <- tax_glom(ps_merged, "Source", NArm = FALSE)
ps0 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))

#Create melted dataframe
df <- psmelt(ps0)

df_FVT <- subset(df, Source %in% "Inoculum")


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

tax <- data.frame(Group=rep("Viruse",length(taxonomy$Source)),
                  Source=taxonomy$Source)

rownames(tax) <- rownames(taxonomy)




# Convert the modified data frame back to a Taxonomy Table
tax <- tax_table(as.matrix(tax))
# Update the original phyloseq object with the new taxonomy table

tax_table(ps_merged) <- tax

ps_merged <- subset_samples(ps_merged, Sample_origin %in% "Pig")

vir.phyl <- tax_glom(ps_merged, "Source", NArm = FALSE)
ps0 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))

#Create melted dataframe
df <- psmelt(ps0)

df_CVT <- subset(df, Source %in% "Inoculum")

##################################1. CVT_MO############################################

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

tax <- data.frame(Group=rep("Viruse",length(taxonomy$Source)),
                  Source=taxonomy$Source)

rownames(tax) <- rownames(taxonomy)




# Convert the modified data frame back to a Taxonomy Table
tax <- tax_table(as.matrix(tax))
# Update the original phyloseq object with the new taxonomy table

tax_table(ps_merged) <- tax

ps_merged <- subset_samples(ps_merged, Sample_origin %in% "Pig")

vir.phyl <- tax_glom(ps_merged, "Source", NArm = FALSE)
ps0 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))

#Create melted dataframe
df <- psmelt(ps0)

df_CVT_MO <- subset(df, Source %in% "Inoculum")

##################################1. CON############################################

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

tax <- data.frame(Group=rep("Viruse",length(taxonomy$Source)),
                  Source=taxonomy$Source)

rownames(tax) <- rownames(taxonomy)




# Convert the modified data frame back to a Taxonomy Table
tax <- tax_table(as.matrix(tax))
# Update the original phyloseq object with the new taxonomy table

tax_table(ps_merged) <- tax

ps_merged <- subset_samples(ps_merged, Sample_origin %in% "Pig")

vir.phyl <- tax_glom(ps_merged, "Source", NArm = FALSE)
ps0 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))

#Create melted dataframe
df <- psmelt(ps0)

df_CON <- subset(df, Source %in% "Inoculum")

##########################MERGE############################

df1 <- rbind(df_FVT,df_CVT,df_CVT_MO,df_CON)


n1 <-  count(df1$sample_Group == "FVT")[[2,2]]
n2 <-  count(df1$sample_Group == "CVT")[[2,2]]
n3 <-  count(df1$sample_Group == "CVT_MO")[[2,2]]
n4 <-  count(df1$sample_Group == "CON")[[2,2]]


stat_bac<- df1 %>%
  wilcox_test(Abundance~sample_Group,
              p.adjust.method = "fdr",
              paired = FALSE, 
              alternative = "two.sided",
              detailed = TRUE) 
stat_bac

df1$sample_Group <- factor(df1$sample_Group, levels =S_groups )
p_share_otu <- ggplot(df1, aes(x= sample_Group, y= Abundance)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.2) + 
  geom_boxplot(outlier.shape = NA, aes(fill=sample_Group)) +
  geom_point(position = position_jitter(w=0.1, h=0),size=2) +
  # labs(color=" NEC")+
  scale_fill_manual(values = simone_cols) +
  scale_x_discrete(labels=c(glue("FVT\nN={n1}"),
                            glue("CVT\nN={n2}"),
                            glue("CVT_MO\nN={n3}"),
                            glue("CON\nN={n4}"))
  )+
  labs(x="", y="Inoculum vOTUs exist in piglet - percentage (Group)", title="Virome") +
  stat_pvalue_manual(stat_bac,label = "p.adj.signif", tip.length = 0, size = 6,
                     y.position = c(NA,NA,NA,NA,NA,95))+
  theme_classic() +
  mytheme_beta
p_share_otu

write_xlsx(stat_bac,"Virome/stat_result/Simone/shared_otu.xlsx")
