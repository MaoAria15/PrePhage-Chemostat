

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################Load sequencing data#######################################

##Import files

#Bacteria

source("2_CSS_phyloseq_bac.R")


############Set data categories - and order



PSB.CSS = prune_samples(sample_sums(PSB) >= 2000, PSB.CSS)
PSB = prune_samples(sample_sums(PSB) >= 2000, PSB)


PSB <- subset_samples(PSB, !Sample_ID %in% "PV2 DN5.2")
PSB.CSS <- subset_samples(PSB.CSS, !Sample_ID %in% "PV2 DN5.2")


#
# saveRDS(PSB, file = "data/corr_data/ps_bacterium.rds")

