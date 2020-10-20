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
  
  carni_res_list[[i]] <- list(cosmos_forward, cosmos_backward)
  names(carni_res_list)[i] <- run_index
  i <- i+1
}

# save(carni_res_list, file = "formated_res_list.Rdata")

# x <- unique(carni_res_list[[1]][[2]]$Nodes)
# y <- unique(carni_res_list[[2]][[2]]$Nodes)

# jaccard_distance(x, y)

jc_vec <- lapply(carni_res_list,function(x,carni_res_list){
  jc <- lapply(carni_res_list, function(y, x){
    x[[1]][[2]]$AvgAct <- as.numeric(x[[1]][[2]]$AvgAct)
    y[[1]][[2]]$AvgAct <- as.numeric(y[[1]][[2]]$AvgAct)
    nodes_A <- unique(x[[1]][[2]][abs(x[[1]][[2]]$AvgAct) > 90 ,"Nodes"])
    nodes_B <- unique(y[[1]][[2]][abs(x[[1]][[2]]$AvgAct) > 90 ,"Nodes"])
    return(jaccard_distance(nodes_A, nodes_B))
  }, x = x)
}, carni_res_list = carni_res_list)

jc_matrix <- as.data.frame(do.call(cbind,jc_vec))
jc_matrix <- jc_matrix[order(as.numeric(row.names(jc_matrix))),]
jc_matrix <- jc_matrix[,order(as.numeric(names(jc_matrix)))]

jc_matrix <- as.data.frame(apply(jc_matrix,2,as.numeric))
row.names(jc_matrix) <- names(jc_matrix)

pheatmap(jc_matrix, cluster_rows = F, cluster_cols = F, display_numbers = T, main = "Forward > 90")

jc_vec <- lapply(carni_res_list,function(x,carni_res_list){
  jc <- lapply(carni_res_list, function(y, x){
    x[[2]][[2]]$AvgAct <- as.numeric(x[[2]][[2]]$AvgAct)
    y[[2]][[2]]$AvgAct <- as.numeric(y[[2]][[2]]$AvgAct)
    nodes_A <- unique(x[[2]][[2]][abs(x[[2]][[2]]$AvgAct) > 90 ,"Nodes"])
    nodes_B <- unique(y[[2]][[2]][abs(x[[2]][[2]]$AvgAct) > 90 ,"Nodes"])
    return(jaccard_distance(nodes_A, nodes_B))
  }, x = x)
}, carni_res_list = carni_res_list)

jc_matrix <- as.data.frame(do.call(cbind,jc_vec))
jc_matrix <- jc_matrix[order(as.numeric(row.names(jc_matrix))),]
jc_matrix <- jc_matrix[,order(as.numeric(names(jc_matrix)))]

jc_matrix <- as.data.frame(apply(jc_matrix,2,as.numeric))
row.names(jc_matrix) <- names(jc_matrix)

pheatmap(jc_matrix, cluster_rows = F, cluster_cols = F, display_numbers = T, main = "Backward > 90")

################

intersect_list <- list()
union_list <- list()
n_cumul_intersect <- c(rep(NA,length(carni_res_list)))
n_cumul_union <- c(rep(NA,length(carni_res_list)))

j <- 1
for(i in as.character(sort(as.numeric(names(carni_res_list)))))
{
  if(i == '0')
  {
    intersect_list[[j]] <- carni_res_list[[i]][[1]][[2]]$Nodes
    union_list[[j]] <- carni_res_list[[i]][[1]][[2]]$Nodes
    n_cumul_intersect[j] <- length(intersect_list[[j]])
    n_cumul_union[j] <- length(intersect_list[[j]])
  } else
  {
    intersect_list[[j]] <- intersect(intersect_list[[j-1]],carni_res_list[[i]][[1]][[2]]$Nodes)
    union_list[[j]] <- union(union_list[[j-1]],carni_res_list[[i]][[1]][[2]]$Nodes)
    n_cumul_intersect[j] <- length(intersect_list[[j]])/length(intersect_list[[1]])
    n_cumul_union[j] <- length(union_list[[j]])/length(union_list[[1]])
  }
  j <- j + 1
}

n_cumul_intersect[1] <- 1
n_cumul_union[1] <- 1

plot(n_cumul_intersect, ylim = c(0,1))

sum(intersect_list[length(intersect_list)][[1]] %in% carni_res_list$`0`[[1]][[2]][carni_res_list$`0`[[1]][[2]]$measured == 1, "Nodes"]) / length(intersect_list[length(intersect_list)][[1]])
