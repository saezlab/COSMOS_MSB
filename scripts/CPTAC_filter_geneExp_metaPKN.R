library(readr)
library(biomaRt)

setwd("~/Dropbox/COSMOS_MSB/")

meta_network_carnival_ready_exch_solved_fullomni <- as.data.frame(
  read_csv("support/metaPKN.csv"))

ttop_CPTAC_CCRCC <- as.data.frame(
  read_csv("data/CPTAC_RNA_ttop.csv"))

expressed_gene_list <- gsub("X","",ttop_CPTAC_CCRCC$ID)

is_expressed <- function(x)
{
  if(!grepl("Metab",x))
  {
    if(gsub("X","",x) %in% expressed_gene_list)
    {
      return(x)
    } else
    {
      if(grepl("XGene[0-9]+__[0-9_]+$",x))
      {
        genes <- gsub("XGene[0-9]+__","",x)
        genes <- strsplit(genes,"_")[[1]]
        if(sum(genes %in% expressed_gene_list) != length(genes))
        {
          return(NA)
        } else
        {
          return(x)
        }
      } else
      {
        if(grepl("XGene[0-9]+__[0-9_]+reverse",x))
        {
          genes <- gsub("XGene[0-9]+__","",x)
          genes <- gsub("_reverse","",genes)
          genes <- strsplit(genes,"_")[[1]]
          if(sum(genes %in% expressed_gene_list) != length(genes))
          {
            return(NA)
          } else
          {
            return(x)
          }
        } else
        {
          return(NA) #changed from NA to x to include reactions not associated with genes. Thanks Vitor.
        }
      }
    } 
  } else
  {
    return(x)
  }
}


# is_expressed("XGene3004__124975_91227")

meta_network_carnival_ready_exch_solved_fullomni$source <- sapply(meta_network_carnival_ready_exch_solved_fullomni$source,is_expressed)
meta_network_carnival_ready_exch_solved_fullomni <- meta_network_carnival_ready_exch_solved_fullomni[complete.cases(meta_network_carnival_ready_exch_solved_fullomni),]
meta_network_carnival_ready_exch_solved_fullomni$target <- sapply(meta_network_carnival_ready_exch_solved_fullomni$target,is_expressed)
meta_network_carnival_ready_exch_solved_fullomni <- meta_network_carnival_ready_exch_solved_fullomni[complete.cases(meta_network_carnival_ready_exch_solved_fullomni),]

write_csv(meta_network_carnival_ready_exch_solved_fullomni,"support/metaPKN_CPTAC_filtered.csv")
