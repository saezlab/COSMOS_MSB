library(readr)

setwd("~/Dropbox/COSMOS_MSB/")

meta_network_carnival_ready_exch_solved_fullomni <- as.data.frame(
  read_csv("support/metaPKN_filtered.csv"))

expressed_gene_list <- as.data.frame(
  read_csv("support/expressed_gene_list.txt", 
           col_names = FALSE))

expressed_gene_list <- as.character(as.vector(as.matrix(expressed_gene_list)))

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
          return(NA)
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

swisscheese_PKN_list <- lapply(seq(2,50,2),function(x,meta_PKN)
  {
  number_missing_edges <- x * length(meta_PKN[,1]) / 100
  holes_to_dig <- sample(1:length(meta_PKN[,1]),number_missing_edges,replace = F)
  swisschesse_PKN <- meta_PKN[-holes_to_dig,]
  
  swisschesse_PKN$source <- sapply(swisschesse_PKN$source,is_expressed)
  swisschesse_PKN <- swisschesse_PKN[complete.cases(swisschesse_PKN),]
  swisschesse_PKN$target <- sapply(swisschesse_PKN$target,is_expressed)
  swisschesse_PKN <- swisschesse_PKN[complete.cases(swisschesse_PKN),]
  
  return(swisschesse_PKN)
}, meta_PKN = meta_network_carnival_ready_exch_solved_fullomni)

names(swisscheese_PKN_list) <- seq(2,50,2)

# save(swisscheese_PKN_list,file = "support/revisions_swisscheese_PKN_list.RData") #The sampling is random. Save your res somewhere safe before overwritting it.
