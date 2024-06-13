setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################Load sequencing data#######################################

##Import files

#Bacteria

source("3_Analysis0_file_loading_and_prep.R")


PSB <- subset_samples(PSB, Group %in% c("FVT","CVT","CVT_MO","CON"))
PSB.CSS <- subset_samples(PSB.CSS, Group %in% c("FVT","CVT","CVT_MO","CON"))