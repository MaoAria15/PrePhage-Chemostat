
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
source("setup.R")

#OTU table
host_mat <- read.csv('Data/Virome/host_pred_formatted.csv', sep = ';',row.names = 1)

otu <- read.delim('Data/Virome/otu_table.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)


otu_classified <- intersect(rownames(otu),rownames(host_mat))
otu_all <- rownames(otu)
missing_votus <- otu_all[!otu_all %in% otu_classified]
Unknown_tax <- data.frame(row.names = missing_votus,
                          Domian = rep("Unknown", length(missing_votus)),
                          Phylum = rep("Unknown", length(missing_votus)),
                          Class = rep("Unknown", length(missing_votus)),
                          Order = rep("Unknown", length(missing_votus)),
                          Family = rep("Unknown", length(missing_votus)),
                          Genus = rep("Unknown", length(missing_votus)),
                          Species = rep("Unknown", length(missing_votus))
                          )

taxonomy <- rbind(Unknown_tax,host_mat)

tax.clean <- data.frame(row.names = row.names(taxonomy),
                        Domian = str_replace(taxonomy[,1], "k__",""),
                        Phylum = str_replace(taxonomy[,2], "p__",""),
                        Class = str_replace(taxonomy[,3], "c__",""),
                        Order = str_replace(taxonomy[,4], "o__",""),
                        Family = str_replace(taxonomy[,5], "f__",""),
                        Genus = str_replace(taxonomy[,6], "g__",""),
                        Species = str_replace(taxonomy[,7], "s__",""),
                        stringsAsFactors = FALSE)

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,7] != ""){
    # tax.clean$Species[i] <- paste(tax.clean$Genus[i], tax.clean$Species[i], sep = " ")
  } else if (tax.clean[i,2] == ""){
    kingdom <- paste("k", tax.clean[i,1], sep = "_")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("p", tax.clean[i,2], sep = "_")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("c", tax.clean[i,3], sep = "_")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("o", tax.clean[i,4], sep = "_")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("f", tax.clean[i,5], sep = "_")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("g ",tax.clean$Genus[i], sep = "_")
  }
}





write.csv(tax.clean,"Data/Virome/host_taxonomy.csv")

