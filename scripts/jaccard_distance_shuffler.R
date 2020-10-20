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

setwd("~/Dropbox/COSMOS_MSB/results/shuffler/")
source("../../scripts/format_cosmos_res_function.R")

jaccard_distance <- function(x, y){length(intersect(x, y))/length(union(x, y))}

files <- list.files(".", recursive = T)
files <- files[grepl("doublerun_res",files)]

ori_file <- "../COSMOS_result/COSMOS_res_session.RData"
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

save(carni_res_list, file = "formated_res_list.Rdata")

x <- unique(carni_res_list[[1]][[2]]$Nodes)
y <- unique(carni_res_list[[2]][[2]]$Nodes)

jaccard_distance(x, y)

jc_vec <- lapply(carni_res_list,function(x,carni_res_list){
  jc <- lapply(carni_res_list, function(y, x){
    nodes_A <- unique(x[[2]]$Nodes)
    nodes_B <- unique(y[[2]]$Nodes)
    return(jaccard_distance(nodes_A, nodes_B))
  }, x = x)
}, carni_res_list = carni_res_list)

jc_matrix <- as.data.frame(do.call(cbind,jc_vec))
jc_matrix <- jc_matrix[order(as.numeric(row.names(jc_matrix))),]
jc_matrix <- jc_matrix[,order(as.numeric(names(jc_matrix)))]

jc_matrix <- as.data.frame(apply(jc_matrix,2,as.numeric))
row.names(jc_matrix) <- names(jc_matrix)

pheatmap(jc_matrix, cluster_rows = F, cluster_cols = F, display_numbers = T)
