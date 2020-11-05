library(readr)
library('org.Hs.eg.db')

setwd("~/Dropbox/COSMOS_MSB/results/")

load("shuffler/2/shuffler_2_carni_doublerun_res_fullomni_met_expfiltered_newDoro_long.RData")

edge_df_frequency <- read_csv("variability/edge_df_frequency.csv")

sif_2p <- readRDS('shuffler/2/full_sif_2p.Rds') 

sif_2p$edgeIDs <- paste0(sif_2p$Node1, '______', sif_2p$Sign, '______', sif_2p$Node2)

# edge_df_frequency <- merge(edge_df_frequency, sif_2p[,c(5,4)], by = 'edgeIDs')

prots <- unique(c(meta_network$source, meta_network$target))
prots <- prots[!(grepl("Metab",prots))]
prots <- gsub("Gene[0-9]+__","",prots)
prots <- gsub("_reverse","",prots)
prots <- gsub("EXCHANGE.*","",prots)
prots <- unique(prots)
prots <- unlist(sapply(prots, function(x){unlist(strsplit(x,split = "_"))}))
entrezid <- gsub('X','',prots)
gene_mapping <- mapIds(org.Hs.eg.db, entrezid, 'SYMBOL', 'ENTREZID')
gene_mapping <- unlist(gene_mapping)
gene_mapping <- gene_mapping[!is.na(gene_mapping)]

meta_network$source <- sapply(meta_network$source, function(x, gene_mapping){
  x <- gsub('^X','',x)
  if(x %in% names(gene_mapping))
  {
    return(gene_mapping[x])
  } else
  {
    return(paste0('X',x))
  }
}, gene_mapping = gene_mapping)

meta_network$target <- sapply(meta_network$target, function(x, gene_mapping){
  x <- gsub('^X','',x)
  if(x %in% names(gene_mapping))
  {
    return(gene_mapping[x])
  } else
  {
    return(paste0('X',x))
  }
}, gene_mapping = gene_mapping)

meta_network$edgeIDs <- paste0(meta_network$source, '______', meta_network$interaction, '______', meta_network$target)
meta_network$ori <- 'present'

edge_df_frequency_in_2pPKN <- merge(edge_df_frequency, meta_network[,c(4,5),drop = F], by = 'edgeIDs', all.x = T)

