setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################Load sequencing data#######################################
##Import files
source("4_simone.R")

level <-"Species"
################Litter 123##############################################################
ps <- subset_samples(PSB, Litter %in% c("1","2","3"))
ps <- tax_glom(ps, level, NArm = FALSE) #select a level to compare


### loop for all the grouping compared with Saline
Groups <- unique(sample_data(ps)$Group)

Groups

Groups <- Groups[Groups!="CVT_MO"]
path_table <- "Bacteriome/stat_result/Simone/Deseq/Litter123/"
dir.create(path_table)

for (group in Groups){
  ps.sub <- prune_samples(sample_data(ps)$Group %in% c(group, "CVT_MO"), ps)
  ps.sub
  # remove all error taxa
  ps.ds <- phyloseq_to_deseq2(ps.sub, ~Group)
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
  write.table(data.frame(sigtab), paste0(path_table,str_replace(group," ",""),"_CVT_MO.tsv"), sep="\t", col.names = NA)
}

################Litter 45##############################################################
ps <- subset_samples(PSB, Litter %in% c("4","5"))
ps <- tax_glom(ps, level, NArm = FALSE) #select a level to compare


### loop for all the grouping compared with Saline
Groups <- unique(sample_data(ps)$Group)

Groups
Groups <- Groups[Groups!="CON"]
path_table <- "Bacteriome/stat_result/Simone/Deseq/Litter45/"
dir.create(path_table)

for (group in Groups){
  ps.sub <- prune_samples(sample_data(ps)$Group %in% c(group, "CON"), ps)
  ps.sub
  # remove all error taxa
  ps.ds <- phyloseq_to_deseq2(ps.sub, ~Group)
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
  write.table(data.frame(sigtab), paste0(path_table,str_replace(group," ",""),"_CON.tsv"), sep="\t", col.names = NA)
}
