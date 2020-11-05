library(readr)

setwd("~/Dropbox/COSMOS_MSB/")

meta_network_carnival_ready_exch_solved_fullomni <- as.data.frame(
  read_csv("support/metaPKN.csv"))

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
# 
# meta_network_carnival_ready_exch_solved_fullomni$source <- sapply(meta_network_carnival_ready_exch_solved_fullomni$source,is_expressed)
# meta_network_carnival_ready_exch_solved_fullomni <- meta_network_carnival_ready_exch_solved_fullomni[complete.cases(meta_network_carnival_ready_exch_solved_fullomni),]
# meta_network_carnival_ready_exch_solved_fullomni$target <- sapply(meta_network_carnival_ready_exch_solved_fullomni$target,is_expressed)
# meta_network_carnival_ready_exch_solved_fullomni <- meta_network_carnival_ready_exch_solved_fullomni[complete.cases(meta_network_carnival_ready_exch_solved_fullomni),]

interaction_suffler_PKN_list <- lapply(rep(2,1),function(x,meta_PKN) #seq(2,50,2)
{
  number_shuffled_edges <- x * length(meta_PKN[,1]) / 100
  
  meta_PKN <- meta_PKN[sample(1:length(meta_PKN[,1]),size = length(meta_PKN[,1]),replace = F),]
  
  interaction_suffler_PKN <- meta_PKN
  
  interaction_suffler_PKN[1:number_shuffled_edges,1] <- interaction_suffler_PKN[sample(1:number_shuffled_edges,size = number_shuffled_edges, replace = F),1]
  
  interaction_suffler_PKN$source <- sapply(interaction_suffler_PKN$source,is_expressed)
  interaction_suffler_PKN <- interaction_suffler_PKN[complete.cases(interaction_suffler_PKN),]
  interaction_suffler_PKN$target <- sapply(interaction_suffler_PKN$target,is_expressed)
  interaction_suffler_PKN <- interaction_suffler_PKN[complete.cases(interaction_suffler_PKN),]
  
  return(interaction_suffler_PKN)
}, meta_PKN = meta_network_carnival_ready_exch_solved_fullomni)

names(interaction_suffler_PKN_list) <- seq(1,1,1)

save(interaction_suffler_PKN_list,file = "support/revisions_interaction_suffler_PKN_2psingle.RData") #The sampling is random. Save your res somewhere safe before overwritting it.

for(i in 1:25)
{
  test  <- interaction_suffler_PKN_list[[i]]
  
  test2 <- meta_network_carnival_ready_exch_solved_fullomni
  
  test <- test[order(test$source),]
  test2 <- test2[order(test2$source),]
  
  test$edge <- paste0(test$source, test$interaction, test$target)
  test2$edge <- paste0(test2$source, test2$interaction, test2$target)
  
  print(sum(test$edge %in% test2$edge))
}
