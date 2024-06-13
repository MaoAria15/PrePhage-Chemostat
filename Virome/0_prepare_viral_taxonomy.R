
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
source("setup.R")

#OTU table
otu_mat <- as.data.frame(read.csv('Data/Virome/taxonomy_lineage.txt',sep = '\t'))

otu <- read.delim('Data/Virome/otu_table.tsv', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)

otu_classified <- intersect(rownames(otu),otu_mat$vOTU)
otu_all <- rownames(otu)
missing_votus <- otu_all[!otu_all %in% otu_classified]
Unknown_tax <- data.frame(vOTU= missing_votus,
                          taxonomy = rep("Unknown;;;;;;;;", length(missing_votus)))

taxonomy <- rbind(Unknown_tax,otu_mat)
rownames <- taxonomy$vOTU
taxonomy <- subset(taxonomy, select = -vOTU)
rownames(taxonomy) <- rownames


tax_clean <- taxonomy %>%
  dplyr::select(taxonomy) %>%
  separate(taxonomy, c("Domain","Realm","Kingdom","Phylum","Class","Order","Family","Subfamily","Genus"), ";")
tax_clean[] <- lapply(tax_clean, function(x) gsub("Unclassified", "", x))

tax_clean <- subset(tax_clean, select = -Subfamily)
tax_clean <- subset(tax_clean, select = -Genus)
tax_clean[] <- lapply(tax_clean, function(x) gsub("Other", "", x))


for (i in 1:nrow(tax_clean)){
  if (tax_clean[i,7] != ""){
    tax_clean$Family[i] <- paste("",tax_clean$Family[i], sep = "")
  }  else if (tax_clean[i,2] == ""){
    domain <- paste("Unclassified", tax_clean[i,1], sep = "_")
    tax_clean[i, 2:7] <- domain
  } else if (tax_clean[i,3] == ""){
    realm <- paste("Unclassified", tax_clean[i,2], sep = "_")
    tax_clean[i, 3:7] <- realm
  } else if (tax_clean[i,4] == ""){
    phylum <- paste("Unclassified", tax_clean[i,3], sep = "_")
    tax_clean[i, 4:7] <- phylum
  } else if (tax_clean[i,5] == ""){
    class <- paste("Unclassified", tax_clean[i,4], sep = "_")
    tax_clean[i, 5:7] <- class
  } else if (tax_clean[i,6] == ""){
    order <- paste("Unclassified", tax_clean[i,5], sep = "_")
    tax_clean[i, 6:7] <- order
  } else if (tax_clean[i,7] == ""){
    tax_clean$Genus[i] <- paste("Unclassified ",tax_clean$Order[i], sep = "_")
  }
}
tax_clean[] <- lapply(tax_clean, function(x) gsub("Unclassified_Unknown", "Unknown", x))







# Still few level need to be maunally fit

write.csv(tax_clean,"Data/Virome/Viral_taxonomy.csv")

