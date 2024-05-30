setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# source("4_Chemostat.R")
source("5_Chemostat_deduct_media.R") #dedect media

ps.rarefied <- subset_samples(PSV.no.realm, Sample_time_point %in% group_point) 

taxonomy <- read.csv('Data/Virome/eukaryotic_taxonomy.csv' ,row.names = 1)
tax_table(ps.rarefied) <- as.matrix(taxonomy)

vir.phyl <- tax_glom(ps.rarefied, "Family", NArm = FALSE)
ps0 <- transform_sample_counts(vir.phyl, function(x) x / sum(x))

#Load phyloseq files to ampvis2 format


PSVamp <- amp_load(ps0)
####################p_timepoint#########################################
PSVamp$metadata$Substrate <- factor(PSVamp$metadata$Substrate, levels = group_subsrate)
PSVamp$metadata$Sample_time_point <- factor(PSVamp$metadata$Sample_time_point, levels = group_point)

p_timepoint <-amp_heatmap(PSVamp,
                                           group_by = "Substrate",
                                           facet_by = "Sample_time_point",
                                           plot_values = FALSE,
                                           tax_show = 12,
                                           tax_aggregate = "Family" ,
                                           tax_add = "Class",
                                           tax_empty = "best",
                                           showRemainingTaxa = TRUE,
                                           normalise = TRUE,
                                           color_vector = c("white", "red"),
                                           plot_colorscale = "sqrt",
                                           plot_legendbreaks = c(1, 10, 25, 50)
) +
  ggtitle(paste("Eukaryotic heatmap")) +
  scale_x_discrete(breaks=group_subsrate,
                   labels=group_subsrate)+
  theme_classic() +
  mytheme_abundance
p_timepoint


#############################p_substrate###########################################


PSVamp$metadata$Substrate <- factor(PSVamp$metadata$Substrate, levels = c("Inoculum","HMO","Lactose"))
PSVamp$metadata$Sample_time_point <- factor(PSVamp$metadata$Sample_time_point, levels = group_point)

p_substrate <-amp_heatmap(PSVamp,
                                  group_by = "Sample_time_point",
                                  facet_by = "Substrate",
                                  plot_values = FALSE,
                                  tax_show = 12,
                                  tax_aggregate = "Family" ,
                                  tax_add = "Class",
                                  tax_empty = "best",
                                  showRemainingTaxa = TRUE,
                                  normalise = TRUE,
                                  color_vector = c("white", "red"),
                                  plot_colorscale = "sqrt",
                                  plot_legendbreaks = c(1, 10, 25, 50)
) +
  ggtitle(paste("Eukaryotic heatmap")) +
  scale_x_discrete(breaks=group_point,
                   labels=group_point)+
  theme_classic() +
  mytheme_abundance
p_substrate
