setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################Load sequencing data#######################################
##Import files
source("4_simone.R")

level <-"Species"
################Litter 123##############################################################
ps <- subset_samples(PSB, Litter %in% c("1","2","3"))
ps <- tax_glom(ps, level, NArm = FALSE) #select a level to compare


### loop for all the grouping compared with Saline
Groups <- unique(sample_data(ps)$NEC_NO2 )

Groups

Groups <- Groups[Groups!="NO"]
path_table <- "Bacteriome/stat_result/Simone/Deseq/Litter123/"
dir.create(path_table)

for (group in Groups){
  ps.sub <- prune_samples(sample_data(ps)$NEC_NO2 %in% c(group, "NO"), ps)
  ps.sub
  # remove all error taxa
  ps.ds <- phyloseq_to_deseq2(ps.sub, ~NEC_NO2)
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
  write.table(data.frame(sigtab), paste0(path_table,str_replace(group," ",""),"_NO_litter123.tsv"), sep="\t", col.names = NA)
}

################Litter 45##############################################################
ps <- subset_samples(PSB, Litter %in% c("4","5"))
ps <- tax_glom(ps, level, NArm = FALSE) #select a level to compare


### loop for all the grouping compared with Saline
Groups <- unique(sample_data(ps)$NEC_NO2 )

Groups

Groups <- Groups[Groups!="NO"]
path_table <- "Bacteriome/stat_result/Simone/Deseq/Litter45_species/"
dir.create(path_table)

for (group in Groups){
  ps.sub <- prune_samples(sample_data(ps)$NEC_NO2 %in% c(group, "NO"), ps)
  ps.sub
  # remove all error taxa
  ps.ds <- phyloseq_to_deseq2(ps.sub, ~NEC_NO2)
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
  write.table(data.frame(sigtab), paste0(path_table,str_replace(group," ",""),"_NO_litter45.tsv"), sep="\t", col.names = NA)
}


################Litter 45##############################################################
ps <- subset_samples(PSB, Litter %in% c("1","2","3","4","5"))
ps <- tax_glom(ps, level, NArm = FALSE) #select a level to compare


### loop for all the grouping compared with Saline
Groups <- unique(sample_data(ps)$NEC_NO2 )

Groups

Groups <- Groups[Groups!="NO"]
path_table <- "Bacteriome/stat_result/Simone/Deseq/"
dir.create(path_table)

for (group in Groups){
  ps.sub <- prune_samples(sample_data(ps)$NEC_NO2 %in% c(group, "NO"), ps)
  ps.sub
  # remove all error taxa
  ps.ds <- phyloseq_to_deseq2(ps.sub, ~NEC_NO2)
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
  write.table(data.frame(sigtab), paste0(path_table,str_replace(group," ",""),"_NO_litter12345.tsv"), sep="\t", col.names = NA)
}
