setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################Load sequencing data#######################################

##Import files

#Virome

source("3_CSS_phyloseq_vir.R")


PSV.no.realm <- subset_samples(PSV.no.realm, Sample_origin %in% c("Chemostat"))
PSV.no.Realm.CSS <- subset_samples(PSV.no.Realm.CSS, Sample_origin %in% c("Chemostat"))
PSV.HOST <- subset_samples(PSV.HOST, Sample_origin %in% c("Chemostat"))

PSV.no.realm <- subset_samples(PSV.no.realm, !Sample_id %in% "PV2 DN5.2")
PSV.no.Realm.CSS <- subset_samples(PSV.no.Realm.CSS, !Sample_id %in% "PV2 DN5.2")
PSV.HOST <- subset_samples(PSV.HOST, !Sample_id %in% "PV2 DN5.2")