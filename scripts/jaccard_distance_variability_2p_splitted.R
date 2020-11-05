####Running simulations as a job array instead of individual jobs seems to impact variability. This results cannot be compared with the rest

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
library(visNetwork)

setwd("~/Dropbox/COSMOS_MSB/results/variability_2p/")
source("../../scripts/format_cosmos_res_function.R")
source("../../scripts/COSMOS_functions.R")
source("../../scripts/piano_function.R")

jaccard_distance <- function(x, y){length(intersect(x, y))/length(union(x, y))}

files <- list.files(".", recursive = T)
files <- files[grepl("doublerun_res",files)]

# ori_file <- "../COSMOS_result/COSMOS_res_session.RData"
# files <- c(ori_file,files)

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

####FORWARD

#NODES

jc_vec <- lapply(carni_res_list,function(x,carni_res_list){
  jc <- lapply(carni_res_list, function(y, x){
    x[[1]][[2]]$AvgAct <- as.numeric(x[[1]][[2]]$AvgAct)
    y[[1]][[2]]$AvgAct <- as.numeric(y[[1]][[2]]$AvgAct)
    nodes_A <- unique(x[[1]][[2]][abs(x[[1]][[2]]$AvgAct) > 0 & x[[1]][[2]]$measured == 0,"Nodes"])
    # nodes_A <- unique(x[[1]][[2]][abs(x[[1]][[2]]$AvgAct) > 0,"Nodes"])
    # measured_A <- unique(x[[1]][[2]][x[[1]][[2]]$measured == 1 ,"Nodes"])
    # print(length(nodes_A))
    # print(length(measured_A))
    nodes_B <- unique(y[[1]][[2]][abs(y[[1]][[2]]$AvgAct) > 0 & y[[1]][[2]]$measured == 0,"Nodes"])
    # nodes_B <- unique(y[[1]][[2]][abs(x[[1]][[2]]$AvgAct) > 0,"Nodes"])
    # measured_B <- unique(y[[1]][[2]][y[[1]][[2]]$measured == 1 ,"Nodes"])
    # print(length(nodes_B))
    # print(length(measured_B))
    return(jaccard_distance(nodes_A, nodes_B))
  }, x = x)
}, carni_res_list = carni_res_list)

jc_matrix <- as.data.frame(do.call(cbind,jc_vec))
jc_matrix <- jc_matrix[order(as.numeric(row.names(jc_matrix))),]
jc_matrix <- jc_matrix[,order(as.numeric(names(jc_matrix)))]

jc_matrix <- as.data.frame(apply(jc_matrix,2,as.numeric))
row.names(jc_matrix) <- names(jc_matrix)
jc_matrix[jc_matrix >= 0.99999] <- NA
mean(unlist(c(jc_matrix)), na.rm = T)

jc_vec <- unlist(c(jc_matrix))
jc_vec <- jc_vec[!is.na(jc_vec)]

plot(density(jc_vec))

pheatmap(jc_matrix, cluster_rows = T, cluster_cols = T, display_numbers = T, main = "Nodes forward")

#EDGES

jc_vec <- lapply(carni_res_list,function(x,carni_res_list){
  jc <- lapply(carni_res_list, function(y, x){
    df_A <- x[[1]][[1]]
    df_A$edgeID <- paste0(x[[1]][[1]][,1],x[[1]][[1]][,2],x[[1]][[1]][,3])
    df_B <- y[[1]][[1]]
    df_B$edgeID <- paste0(y[[1]][[1]][,1],y[[1]][[1]][,2],y[[1]][[1]][,3])
    return(jaccard_distance(df_A$edgeID, df_B$edgeID))
  }, x = x)
}, carni_res_list = carni_res_list)

jc_matrix <- as.data.frame(do.call(cbind,jc_vec))
jc_matrix <- jc_matrix[order(as.numeric(row.names(jc_matrix))),]
jc_matrix <- jc_matrix[,order(as.numeric(names(jc_matrix)))]

jc_matrix <- as.data.frame(apply(jc_matrix,2,as.numeric))
row.names(jc_matrix) <- names(jc_matrix)
jc_matrix[jc_matrix >= 0.99999] <- NA
mean(unlist(c(jc_matrix)), na.rm = T)

jc_vec <- unlist(c(jc_matrix))
jc_vec <- jc_vec[!is.na(jc_vec)]

plot(density(jc_vec))

pheatmap(jc_matrix, cluster_rows = T, cluster_cols = T, display_numbers = T, main = "Edges forward")

####BACKWARD

#NODES
jc_vec <- lapply(carni_res_list,function(x,carni_res_list){
  jc <- lapply(carni_res_list, function(y, x){
    x[[2]][[2]]$AvgAct <- as.numeric(x[[2]][[2]]$AvgAct)
    y[[2]][[2]]$AvgAct <- as.numeric(y[[2]][[2]]$AvgAct)
    nodes_A <- unique(x[[2]][[2]][abs(x[[2]][[2]]$AvgAct) > 0 & x[[2]][[2]]$measured == 0,"Nodes"])
    # nodes_A <- unique(x[[2]][[2]][abs(x[[2]][[2]]$AvgAct) > 0,"Nodes"])
    # measured_A <- unique(x[[2]][[2]][x[[2]][[2]]$measured == 1 ,"Nodes"])
    # print(length(nodes_A))
    # print(length(measured_A))
    nodes_B <- unique(y[[2]][[2]][abs(y[[2]][[2]]$AvgAct) > 0 & y[[2]][[2]]$measured == 0,"Nodes"])
    # nodes_B <- unique(y[[2]][[2]][abs(x[[2]][[2]]$AvgAct) > 0,"Nodes"])
    # measured_B <- unique(y[[2]][[2]][y[[2]][[2]]$measured == 1 ,"Nodes"])
    return(jaccard_distance(nodes_A, nodes_B))
  }, x = x)
}, carni_res_list = carni_res_list)

jc_matrix <- as.data.frame(do.call(cbind,jc_vec))
jc_matrix <- jc_matrix[order(as.numeric(row.names(jc_matrix))),]
jc_matrix <- jc_matrix[,order(as.numeric(names(jc_matrix)))]

jc_matrix <- as.data.frame(apply(jc_matrix,2,as.numeric))
row.names(jc_matrix) <- names(jc_matrix)
jc_matrix[jc_matrix >= 0.999] <- NA
mean(unlist(c(jc_matrix)), na.rm = T)
pheatmap(jc_matrix, cluster_rows = T, cluster_cols = T, display_numbers = T, main = "Nodes backward")

#EDGES

jc_vec <- lapply(carni_res_list,function(x,carni_res_list){
  jc <- lapply(carni_res_list, function(y, x){
    df_A <- x[[2]][[1]]
    df_A$edgeID <- paste0(x[[2]][[1]][,1],x[[2]][[1]][,2],x[[2]][[1]][,3])
    df_B <- y[[2]][[1]]
    df_B$edgeID <- paste0(y[[2]][[1]][,1],y[[2]][[1]][,2],y[[2]][[1]][,3])
    return(jaccard_distance(df_A$edgeID, df_B$edgeID))
  }, x = x)
}, carni_res_list = carni_res_list)

jc_matrix <- as.data.frame(do.call(cbind,jc_vec))
jc_matrix <- jc_matrix[order(as.numeric(row.names(jc_matrix))),]
jc_matrix <- jc_matrix[,order(as.numeric(names(jc_matrix)))]

jc_matrix <- as.data.frame(apply(jc_matrix,2,as.numeric))
row.names(jc_matrix) <- names(jc_matrix)
jc_matrix[jc_matrix >= 0.999] <- NA
mean(unlist(c(jc_matrix)), na.rm = T)

pheatmap(jc_matrix, cluster_rows = T, cluster_cols = T, display_numbers = T, main = "Edge backward")

################

##With input nodes
#Forward
intersect_list <- list()
union_list <- list()
n_cumul_intersect <- c(rep(NA,length(carni_res_list)))
n_cumul_union <- c(rep(NA,length(carni_res_list)))


j <- 1
for(i in as.character(sort(as.numeric(names(carni_res_list)))))
{
  if(j == 1)
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

#Backward
intersect_list <- list()
union_list <- list()
n_cumul_intersect <- c(rep(NA,length(carni_res_list)))
n_cumul_union <- c(rep(NA,length(carni_res_list)))


j <- 1
for(i in as.character(sort(as.numeric(names(carni_res_list)))))
{
  if(j == 1)
  {
    intersect_list[[j]] <- carni_res_list[[i]][[2]][[2]]$Nodes
    union_list[[j]] <- carni_res_list[[i]][[2]][[2]]$Nodes
    n_cumul_intersect[j] <- length(intersect_list[[j]])
    n_cumul_union[j] <- length(intersect_list[[j]])
  } else
  {
    intersect_list[[j]] <- intersect(intersect_list[[j-1]],carni_res_list[[i]][[2]][[2]]$Nodes)
    union_list[[j]] <- union(union_list[[j-1]],carni_res_list[[i]][[2]][[2]]$Nodes)
    n_cumul_intersect[j] <- length(intersect_list[[j]])/length(intersect_list[[1]])
    n_cumul_union[j] <- length(union_list[[j]])/length(union_list[[1]])
  }
  j <- j + 1
}

plot(n_cumul_intersect, ylim = c(0,1))
n_cumul_intersect[j-1]

##Without input nodes

#Forward
intersect_list <- list()
union_list <- list()
n_cumul_intersect <- c(rep(NA,length(carni_res_list)))
n_cumul_union <- c(rep(NA,length(carni_res_list)))

j <- 1
for(i in as.character(sort(as.numeric(names(carni_res_list)))))
{
  if(j == 1)
  {
    intersect_list[[j]] <- carni_res_list[[i]][[1]][[2]][abs(as.numeric(carni_res_list[[i]][[1]][[2]]$AvgAct)) > 0 & carni_res_list[[i]][[1]][[2]]$measured == 0,"Nodes"]
    union_list[[j]] <- carni_res_list[[i]][[1]][[2]][abs(as.numeric(carni_res_list[[i]][[1]][[2]]$AvgAct)) > 0 & carni_res_list[[i]][[1]][[2]]$measured == 0,"Nodes"]
    n_cumul_intersect[j] <- length(intersect_list[[j]])
    n_cumul_union[j] <- length(intersect_list[[j]])
  } else
  {
    intersect_list[[j]] <- intersect(intersect_list[[j-1]],carni_res_list[[i]][[1]][[2]][abs(as.numeric(carni_res_list[[i]][[1]][[2]]$AvgAct)) > 0 & carni_res_list[[i]][[1]][[2]]$measured == 0,"Nodes"])
    union_list[[j]] <- union(union_list[[j-1]],carni_res_list[[i]][[1]][[2]][abs(as.numeric(carni_res_list[[i]][[1]][[2]]$AvgAct)) > 0 & carni_res_list[[i]][[1]][[2]]$measured == 0,"Nodes"])
    n_cumul_intersect[j] <- length(intersect_list[[j]])/length(intersect_list[[1]])
    n_cumul_union[j] <- length(union_list[[j]])/length(union_list[[1]])
  }
  j <- j + 1
}
n_cumul_intersect[1] <- 1
n_cumul_union[1] <- 1

plot(n_cumul_intersect, ylim = c(0,1))
n_cumul_intersect[j-1]


#Backward
intersect_list <- list()
union_list <- list()
n_cumul_intersect <- c(rep(NA,length(carni_res_list)))
n_cumul_union <- c(rep(NA,length(carni_res_list)))

j <- 1
for(i in as.character(sort(as.numeric(names(carni_res_list)))))
{
  if(j == 1)
  {
    intersect_list[[j]] <- carni_res_list[[i]][[2]][[2]][abs(as.numeric(carni_res_list[[i]][[2]][[2]]$AvgAct)) > 00 & carni_res_list[[i]][[2]][[2]]$measured == 0,"Nodes"]
    union_list[[j]] <- carni_res_list[[i]][[2]][[2]][abs(as.numeric(carni_res_list[[i]][[2]][[2]]$AvgAct)) > 00 & carni_res_list[[i]][[2]][[2]]$measured == 0,"Nodes"]
    n_cumul_intersect[j] <- length(intersect_list[[j]])
    n_cumul_union[j] <- length(intersect_list[[j]])
  } else
  {
    intersect_list[[j]] <- intersect(intersect_list[[j-1]],carni_res_list[[i]][[2]][[2]][abs(as.numeric(carni_res_list[[i]][[2]][[2]]$AvgAct)) > 00 & carni_res_list[[i]][[2]][[2]]$measured == 0,"Nodes"])
    union_list[[j]] <- union(union_list[[j-1]],carni_res_list[[i]][[2]][[2]][abs(as.numeric(carni_res_list[[i]][[2]][[2]]$AvgAct)) > 00 & carni_res_list[[i]][[2]][[2]]$measured == 0,"Nodes"])
    n_cumul_intersect[j] <- length(intersect_list[[j]])/length(intersect_list[[1]])
    n_cumul_union[j] <- length(union_list[[j]])/length(union_list[[1]])
  }
  j <- j + 1
}
n_cumul_intersect[1] <- 1
n_cumul_union[1] <- 1
n_cumul_intersect[j-1]

plot(n_cumul_intersect, ylim = c(0,1))

##########################
#COMBINE NETWORK
edge_lists_for <- list() 
edge_lists_back <- list()
for(i in 1:length(carni_res_list))
{
  if(i == 1)
  {
    intersect_network_forward <- carni_res_list[[i]][[1]][[1]]
    intersect_network_forward$edgeID <- paste0(intersect_network_forward[,1],'______',intersect_network_forward[,2],'______',intersect_network_forward[,3])
    edge_lists_for[[i]] <- intersect_network_forward$edgeID
    
    intersect_network_backward <- carni_res_list[[i]][[2]][[1]]
    intersect_network_backward$edgeID <- paste0(intersect_network_backward[,1],'______',intersect_network_backward[,2],'______',intersect_network_backward[,3])
    edge_lists_back[[i]] <- intersect_network_backward$edgeID
  } else
  {
    new_net_for <- carni_res_list[[i]][[1]][[1]]
    new_net_for$edgeID <- paste0(new_net_for[,1],'______',new_net_for[,2],'______',new_net_for[,3])
    edge_lists_for[[i]] <- new_net_for$edgeID
    intersect_network_forward <- intersect_network_forward[intersect_network_forward$edgeID %in% new_net_for$edgeID,]
    
    new_net_back <- carni_res_list[[i]][[2]][[1]]
    new_net_back$edgeID <- paste0(new_net_back[,1],'______',new_net_back[,2],'______',new_net_back[,3])
    edge_lists_back[[i]] <- new_net_back$edgeID
    intersect_network_backward <- intersect_network_backward[intersect_network_backward$edgeID %in% new_net_back$edgeID,]
  }
  # print(length(intersect_network_forward[,1]))
  # print(length(intersect_network_backward[,1]))
  # print('STOP')
}

#edge profiles forward
edge_lists_for_shuffled <- lapply(1:100, function(x,edge_lists_for)
{
  return(edge_lists_for[sample(1:length(edge_lists_for), length(edge_lists_for), replace = F)])
}, edge_lists_for = edge_lists_for)

edge_lists_for <- 
  
  edge_df_for <- as.data.frame(matrix(NA, length(unique(unlist(edge_lists_for))), length(edge_lists_for)))
edge_df_for$edgeID <- unique(unlist(edge_lists_for))

for(i in 1:length(edge_lists_for))
{
  edge_df_for[,i] <- as.character(edge_df_for$edgeID) %in% as.character(edge_lists_for[[i]])
}
row.names(edge_df_for) <- edge_df_for$edgeID
edge_df_for <- edge_df_for[,-length(edge_df_for[1,])]
edge_df_for[edge_df_for == T] <- 100 

pheatmap(edge_df_for)

edge_df_for_cumul <- edge_df_for
edge_df_for_reldif <- edge_df_for
for(i in 1:length(edge_df_for_cumul[1,]))
{
  if(i > 1)
  {
    edge_df_for_cumul[,i] <- as.numeric(rowMeans(edge_df_for[,1:i]))
    edge_df_for_reldif[,i] <- abs(edge_df_for_cumul[,i] - edge_df_for_cumul[,i-1]) / as.numeric(rowMeans(edge_df_for_cumul[,c(i-1,i)]))
    edge_df_for_reldif[is.nan(edge_df_for_reldif[,i]),i] <- 0
  }
}

edge_df_for_reldif[,1] <- 0
average_reldifs_for <- apply(edge_df_for_reldif, 2, mean)

plot(average_reldifs_for[-1])

#edge profiles backward
edge_df_back <- as.data.frame(matrix(NA, length(unique(unlist(edge_lists_back))), length(edge_lists_back)))
edge_df_back$edgeID <- unique(unlist(edge_lists_back))

for(i in 1:length(edge_lists_back))
{
  edge_df_back[,i] <- as.character(edge_df_back$edgeID) %in% as.character(edge_lists_back[[i]])
}
row.names(edge_df_back) <- edge_df_back$edgeID
edge_df_back <- edge_df_back[,-length(edge_df_back[1,])]
edge_df_back[edge_df_back == T] <- 100 

pheatmap(edge_df_back)

edge_df_back_cumul <- edge_df_back
edge_df_back_reldif <- edge_df_back
for(i in 1:length(edge_df_back_cumul[1,]))
{
  if(i > 1)
  {
    edge_df_back_cumul[,i] <- as.numeric(rowMeans(edge_df_back[,1:i]))
    edge_df_back_reldif[,i] <- abs(edge_df_back_cumul[,i] - edge_df_back_cumul[,i-1]) / as.numeric(rowMeans(edge_df_back_cumul[,c(i-1,i)]))
    edge_df_back_reldif[is.nan(edge_df_back_reldif[,i]),i] <- 0
  }
}

edge_df_back_reldif[,1] <- 0
average_reldifs_back <- apply(edge_df_back_reldif, 2, mean)

plot(average_reldifs_back[-1])

######
intersect_att_for <- carni_res_list[[i]][[1]][[2]]
intersect_att_for <- intersect_att_for[intersect_att_for$Nodes %in% intersect_network_forward[,1] | intersect_att_for$Nodes %in% intersect_network_forward[,3],]

intersect_att_back <- carni_res_list[[i]][[2]][[2]]
intersect_att_back <- intersect_att_back[intersect_att_back$Nodes %in% intersect_network_backward[,1] | intersect_att_back$Nodes %in% intersect_network_backward[,3],]

test_for <- intersect_att_for[intersect_att_for$Nodes %in% intersect_att_back$Nodes,]
test_back <- intersect_att_back[intersect_att_back$Nodes %in% intersect_att_for$Nodes,]

test_for$nodeID <- paste0(test_for$Nodes,test_for$Activity)
test_back$nodeID <- paste0(test_back$Nodes, test_back$Activity)

print(test_for[!(test_for$nodeID %in% test_back$nodeID),])
print(test_back[!(test_back$nodeID %in% test_for$nodeID),])

intersect_network_full <- as.data.frame(rbind(intersect_network_forward, intersect_network_backward))
intersect_network_full <- intersect_network_full[!duplicated(intersect_network_full$edgeID),]

intersect_att_full <- as.data.frame(rbind(intersect_att_for, intersect_att_back))
intersect_att_full$nodeID <- paste0(intersect_att_full$Nodes, intersect_att_full$Activity)
intersect_att_full <- intersect_att_full[!duplicated(intersect_att_full$nodeID),]

write_csv(intersect_network_full, 'intersect_network_full_sif.csv')
write_csv(intersect_att_full, 'intersect_att_full.csv')

sif_1 <- carni_res_list[['0']][[1]][[1]]
nodes_1 <- carni_res_list[['0']][[1]][[2]]

display_node_neighboorhood("BCAT1",sif = sif_1[,c(1,3,2,4)], att = nodes_1, n = 3)

sif_1 <- carni_res_list[[4]][[1]][[1]]
nodes_1 <- carni_res_list[[4]][[1]][[2]]

display_node_neighboorhood("BCAT1",sif = sif_1[,c(1,3,2,4)], att = nodes_1, n = 3)


write_csv(intersect_network_full, 'intersect_network_full_sif.csv')
write_csv(intersect_att_full, 'intersect_att_full.csv')