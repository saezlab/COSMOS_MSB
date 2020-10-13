library(CARNIVAL) # load CARNIVAL library
library(biomaRt)
library(gsubfn)
library(readr)
library(rpubchem)
library(stringr)
library(dplyr)
library(dorothea)
library(OmnipathR)
library(igraph)
library(jaccard)
library(pheatmap)

setwd("~/Dropbox/COSMOS_MSB/results/input_range/")

source("../../scripts/format_cosmos_res_function.R")

jaccard_distance <- function(x, y){length(intersect(x, y))/length(union(x, y))}

files <- list.files(".", recursive = T)
files <- files[grepl("doublerun_res",files)]

ori_file <- "../COSMOS_res_session.RData"
files <- c(ori_file,files)

metab_to_pubchem <- as.data.frame(read_csv("../../support/metab_to_pubchem.csv"))



carni_res_list <- list()
i <- 1
for(file in files)
{
  print(i)
  load(file)
  cosmos_forward <- format_COSMOS_res(cosmos_res = CARNIVAL_Result_rerun,
                                      metab_mapping = metab_to_pubchem,
                                      measured_nodes = c(names(metab_input_carnival),names(signaling_input_carnival)),
                                      omnipath_ptm = read_rds("../../support/omnipath_ptm.rds"), gene_mapping = "else")
  cosmos_backward <- format_COSMOS_res(cosmos_res = CARNIVAL_Result_2,
                                       metab_mapping = metab_to_pubchem,
                                       measured_nodes = c(names(metab_input_carnival),names(signaling_input_carnival)),
                                       omnipath_ptm = read_rds("../../support/omnipath_ptm.rds"), gene_mapping = "else")
  run_index <- gsub("[/].*","",file)
  run_index <- gsub("[.][.]","0",run_index)
  
  full_sif <- as.data.frame(unique(rbind(cosmos_forward[[1]],cosmos_backward[[1]])))
  full_att <- as.data.frame(unique(rbind(cosmos_forward[[2]][,-6],cosmos_backward[[2]][,-6])))
  full_att <- unique(full_att[,c(1,6,7,8)])
  full_att[duplicated(full_att$Nodes),]
  
  full_att$Nodes <- gsub(",","_",full_att$Nodes)
  full_sif$Node1 <- gsub(",","_",full_sif$Node1)
  full_sif$Node2 <- gsub(",","_",full_sif$Node2)
  full_sif$Weight <- as.numeric(full_sif$Weight)
  
  full_sif <- full_sif %>% 
    group_by(Node1, Node2, Sign) %>% 
    summarise(Weight= sum(Weight))
  
  full_sif[full_sif$Weight > 100, "Weight"] <- 100
  
  full_att <- full_att[full_att$Nodes %in% full_sif$Node1 | full_att$Nodes %in% full_sif$Node2,]
  
  carni_res_list[[i]] <- list(full_sif, full_att)
  names(carni_res_list)[i] <- run_index
  i <- i+1
}

source("../../scripts/COSMOS_functions.R")
source("../../scripts/piano_function.R")

#Load inputs
meta_network <- as.data.frame(
  read_csv("../../support/metaPKN_filtered.csv"))

prots <- unique(c(meta_network$source,meta_network$target))
prots <- prots[!grepl("XMetab",prots)]
prots <- gsub("^X","",prots)
prots <- gsub("Gene[0-9]+__","",prots)
prots <- gsub("_reverse","",prots)
prots <- gsub("EXCHANGE.*","",prots)
prots <- unique(prots)
prots <- unlist(sapply(prots, function(x){unlist(strsplit(x,split = "_"))}))

gene_mapping <- "else"

if(gene_mapping == "ensembl")
{
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  
  G_list <- getBM(filters = "entrezgene_id",
                  attributes = c('hgnc_symbol','entrezgene_id', "description"),
                  values = prots, mart = ensembl)
  gene_mapping <- G_list[,1]
  names(gene_mapping) <- G_list[,2]
  
} else
{
  entrezid <- prots
  gene_mapping <- mapIds(org.Hs.eg.db, entrezid, 'SYMBOL', 'ENTREZID')
  gene_mapping <- unlist(gene_mapping)
  gene_mapping <- gene_mapping[!is.na(gene_mapping)]
}



background <- sapply(prots, function(x, translation_dictionary){
  return(translation_dictionary[x])
},translation_dictionary = gene_mapping)
names(background) <- 1:length(background)

background <- unique(background)

hallmarks <- gmt_to_df("../../support/h.all.v7.1.symbols.gmt_20200709.gmt")

ORA_res_list <- lapply(carni_res_list, function(x)
{
  nodes <- unique(x[[2]]$Nodes)
  ora_res <- runGSAhyper(nodes, universe = background, gsc = loadGSC(hallmarks))
  return(ora_res)
})

oraCor_vec <- lapply(ORA_res_list,function(x,ORA_res_list){
  jc <- lapply(ORA_res_list, function(y, x){
    ora_A <- as.data.frame(x$pvalues)
    ora_B <- as.data.frame(y$pvalues)
    
    ora_A$pathway <- row.names(ora_A)
    ora_B$pathway <- row.names(ora_B)
    
    ora_comp <- merge(ora_A, ora_B, by = "pathway")
    
    return(cor(ora_comp[,2], ora_comp[,3]))
  }, x = x)
}, ORA_res_list = ORA_res_list)

oraCor_matrix <- as.data.frame(do.call(cbind,oraCor_vec))
oraCor_matrix <- oraCor_matrix[order(as.numeric(row.names(oraCor_matrix))),]
oraCor_matrix <- oraCor_matrix[,order(as.numeric(names(oraCor_matrix)))]

oraCor_matrix <- as.data.frame(apply(oraCor_matrix,2,as.numeric))
row.names(oraCor_matrix) <- names(oraCor_matrix)

pheatmap(oraCor_matrix, cluster_rows = F, cluster_cols = F, display_numbers = T)
