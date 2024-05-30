setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################Load sequencing data#######################################

##Import files

#Virome

source("3_CSS_phyloseq_vir.R")


PSV.no.realm <- subset_samples(PSV.no.realm, Sample_origin %in% c("Chemostat_donor","Pig_inoculum"))
PSV.no.Realm.CSS <- subset_samples(PSV.no.Realm.CSS, Sample_origin %in% c("Chemostat_donor","Pig_inoculum"))
PSV.HOST <- subset_samples(PSV.HOST, Sample_origin %in% c("Chemostat_donor","Pig_inoculum"))
