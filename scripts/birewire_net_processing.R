library(readr)

setwd("~/Dropbox/COSMOS_MSB/")

pkn_list <- readRDS("support/shuffled_dsg_list.Rds")

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

filtered_PKN_list <- lapply(pkn_list,function(x)
{
  PKN <- x
  PKN$source <- sapply(as.character(PKN$source),is_expressed)
  PKN <- PKN[complete.cases(PKN),]
  PKN$target <- sapply(as.character(PKN$target),is_expressed)
  PKN <- PKN[complete.cases(PKN),]
  
  return(PKN)
})

names(filtered_PKN_list) <- seq(2,50,2)

filtered_PKN_list <- lapply(filtered_PKN_list,function(x)
{
  PKN <- x
  PKN$sign <- ifelse(as.character(PKN$sign) == "+", 1, -1)
  names(PKN) <- c("source","interaction","target")
  return(PKN)
})

save(filtered_PKN_list,file = "support/berewire_PKN_list.RData") #The sampling is random. Save your res somewhere safe before overwritting it.

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
