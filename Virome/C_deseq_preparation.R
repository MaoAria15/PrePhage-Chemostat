setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################Load sequencing data#######################################
##Import files
# source("4_Chemostat.R")
source("5_Chemostat_deduct_media.R") #dedect media

################Litter 123##############################################################
ps <- subset_samples(PSV.no.realm, Sample_time_point %in% c("Chemostat4"))


### loop for all the grouping compared with Saline
Groups <- unique(sample_data(ps)$Substrate)

Groups

Groups <- Groups[Groups!="HMO"]
path_table <- "Virome/stat_result/Chemostat/Deseq/"
dir.create(path_table)

for (group in Groups){
  ps.sub <- prune_samples(sample_data(ps)$Substrate %in% c(group, "HMO"), ps)
  ps.sub
  # remove all error taxa
  ps.ds <- phyloseq_to_deseq2(ps.sub, ~Substrate)
  # solve rows without a zero, deseq need to calculate the geometric zero, 
  cts <- counts(ps.ds)
  geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
  dds <- estimateSizeFactors(ps.ds, geoMeans=geoMeans)
  ps.ds <-  DESeq(dds, test="Wald", fitType="parametric")
  # result
  res = results(ps.ds, cooksCutoff = FALSE)
  #alpha = 0.0001
  sigtab = res
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
  head(sigtab)
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  write.table(data.frame(sigtab), paste0(path_table,str_replace(group," ",""),"_HMO_chemostat4(deduct_media).tsv"), sep="\t", col.names = NA)
}

