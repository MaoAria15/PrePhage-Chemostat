setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################Load sequencing data#######################################

##Import files

#Bacteria

source("3_Analysis0_file_loading_and_prep.R")


PSB <- subset_samples(PSB, Sample_origin %in% c("Chemostat"))
PSB.CSS <- subset_samples(PSB.CSS, Sample_origin %in% c("Chemostat"))
