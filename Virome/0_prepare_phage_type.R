setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
source("setup.R")

#OTU table
phatype <- read.csv('Data/Virome/phatyp_prediction.csv', sep = ",")
otu <- read.delim('Data/Virome/otu_table.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)


common_row <- intersect(phatype$Accession, rownames(otu))

rownames(phatype) <- phatype$Accession

phatype <- phatype %>% select(-Accession, -Length, -Score)


# Still few level need to be maunally fit

write.csv(phatype,"Data/Virome//Phage_lifestyle.csv")

