library("DOSE")
library("reticulate")
library("enrichplot")
library("ggplot2")
library("clusterProfiler")
library("org.Hs.eg.db")
library("GOSemSim")
library("AnnotationHub")
np <- import("numpy")

library("DOSE")
library("reticulate")
library("enrichplot")
library("ggplot2")
library("GOplot")

setwd("~/Users/johgeist/Documents/AG_Luck/3did_project/revisions/GO_term_enrichment")

# Load their list of genes, without filtering for any expression values:

data(geneList)
genes <- names(geneList) 
df <- data.frame(entrez_id = gene)
write.csv(df,file='DOSE_genes_list.csv') 

# list the remaining clusters file names
list_files <- list.files(path=getwd(), pattern='CLUSTER.npy')
# cluster names
col_names <- sub(".npy",  "", list_files)
# load the clusters into a list:
clusters <- lapply(list_files, function(x) np$load(x))
# change the names of the list:
names(clusters) <- col_names

# I will do GO-BP, GO-CC, GO-MF and DO enrichment
# the background are all genes in HuRI

BP_enrichment <- lapply(clusters, function(x) enrichGO(gene = x, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = TRUE))
BP_enrichment_dataframes <- lapply(BP_enrichment, function(x) data.frame(ID = x$ID, Description = x$Description, GeneRatio = x$GeneRatio, BgRatio = x$BgRatio, pvalue = x$pvalue, p.adjust = x$p.adjust, qvalue = x$qvalue, geneID = x$geneID, Count = x$Count))
for (i in seq_along(BP_enrichment_dataframes)){
  write.table(BP_enrichment_dataframes[[i]], file=sprintf('BP_%s_general_background.csv', names(BP_enrichment_dataframes)[i]), row.names = FALSE)
}

CC_enrichment <- lapply(clusters, function(x) enrichGO(gene = x, OrgDb = org.Hs.eg.db, ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = TRUE))
CC_enrichment_dataframes <- lapply(CC_enrichment, function(x) data.frame(ID = x$ID, Description = x$Description, GeneRatio = x$GeneRatio, BgRatio = x$BgRatio, pvalue = x$pvalue, p.adjust = x$p.adjust, qvalue = x$qvalue, geneID = x$geneID, Count = x$Count))
for (i in seq_along(CC_enrichment_dataframes)){
  write.table(CC_enrichment_dataframes[[i]], file=sprintf('CC_%s_general_background.csv', names(CC_enrichment_dataframes)[i]), row.names = FALSE)
}

MF_enrichment <- lapply(clusters, function(x) enrichGO(gene = x, OrgDb = org.Hs.eg.db, ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = TRUE))
MF_enrichment_dataframes <- lapply(MF_enrichment, function(x) data.frame(ID = x$ID, Description = x$Description, GeneRatio = x$GeneRatio, BgRatio = x$BgRatio, pvalue = x$pvalue, p.adjust = x$p.adjust, qvalue = x$qvalue, geneID = x$geneID, Count = x$Count))
for (i in seq_along(MF_enrichment_dataframes)){
  write.table(MF_enrichment_dataframes[[i]], file=sprintf('MF_%s_general_background.csv', names(MF_enrichment_dataframes)[i]), row.names = FALSE)
}

DO_enrichment <- lapply(clusters, function(x) enrichDO(x, ont="DO", pvalueCutoff = 0.05, pAdjustMethod="BH", minGSSize=5, readable=TRUE))
DO_enrichment_dataframes <- lapply(DO_enrichment, function(x) data.frame(ID = x$ID, Description = x$Description, GeneRatio = x$GeneRatio, BgRatio = x$BgRatio, pvalue = x$pvalue, p.adjust = x$p.adjust, qvalue = x$qvalue, geneID = x$geneID, Count = x$Count))
for (i in seq_along(DO_enrichment_dataframes)){
  write.table(DO_enrichment_dataframes[[i]], file=sprintf('DO_%s_general_background.csv', names(DO_enrichment_dataframes)[i]), row.names = FALSE)
}

# I will do GO-BP, GO-CC, GO-MF and DO enrichment for every cluster


