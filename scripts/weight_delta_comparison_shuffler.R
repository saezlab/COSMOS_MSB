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
library(ggplot2)
library(reshape2)

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
  
  carni_res_list[[i]] <- list(cosmos_forward, cosmos_backward, run_index)
  names(carni_res_list)[i] <- run_index
  i <- i+1
}

ori_carni_res <- carni_res_list[['0']]

carni_res_list <- carni_res_list[-1]

delta_list_forward <- lapply(carni_res_list, function(x, ori_carni_res){
  forward_shuffle <- x[[1]][[1]]
  forward_shuffle$Weight <- as.numeric(forward_shuffle$Weight)
  forward_shuffle$edgeIDs <- paste0(forward_shuffle$Node1, '______', forward_shuffle$Sign, '______', forward_shuffle$Node2)
  forward_shuffle <- forward_shuffle[,c(5,4)]
  forward_shuffle <- aggregate(Weight ~ edgeIDs, forward_shuffle, sum)
  forward_shuffle$Weight <- ifelse(forward_shuffle$Weight  > 100, 100, forward_shuffle$Weight)
  
  
  forward_ori <- ori_carni_res[[1]][[1]]
  forward_ori$Weight <- as.numeric(forward_ori$Weight)
  forward_ori$edgeIDs <- paste0(forward_ori$Node1, '______', forward_ori$Sign, '______', forward_ori$Node2)
  forward_ori <- forward_ori[,c(5,4)]
  forward_ori <- aggregate(Weight ~ edgeIDs, forward_ori, sum)
  forward_ori$Weight <- ifelse(forward_ori$Weight  > 100, 100, forward_ori$Weight)
  
  delta_df <- merge(forward_ori, forward_shuffle, by = 'edgeIDs', all = T)
  delta_df[is.na(delta_df)] <- 0
  delta_df$delta <- abs(delta_df$Weight.x - delta_df$Weight.y)
  delta_df$index <- x[[3]]
  
  return(delta_df)
}, ori_carni_res = ori_carni_res)

delta_df_for <- as.data.frame(do.call(rbind, delta_list_forward))
delta_df_for$index <- as.numeric(delta_df_for$index)
delta_df_for <- delta_df_for[order(delta_df_for$index, decreasing = F),]

delta_df_for <- delta_df_for[delta_df_for$index %in% c(2,10,20,30,40,50),]

delta_df_for$index <- factor(delta_df_for$index, levels = unique(delta_df_for$index))

median(delta_df_for[delta_df_for$index == "2",'delta'])
p2_totally_different <- delta_df_for[delta_df_for$index == "2" & delta_df_for$delta == 100,]
length(p2_totally_different[,1])/length(delta_df_for[delta_df_for$index == "2",'delta'])

ggplot(delta_df_for, aes(x = index, y = delta, group = index, fill = index)) + 
  geom_boxplot(fill = "lightblue") +
  geom_jitter(alpha = 0.2) +
  stat_summary(fun.y=median, geom="point", shape=23, fill = 'lightgrey', size = 4) +
  theme_minimal() +
  ggtitle('Forward COSMOS weight comparison') +
  xlab('% of edges shuffled') +
  ylab("Absolute edge weight difference")

ggplot(delta_df_for, aes(x = delta, fill = index)) + 
  geom_density(alpha = 0.33) +
  theme_minimal()


#####BACKWARD

delta_list_backward <- lapply(carni_res_list, function(x, ori_carni_res){
  backward_shuffle <- x[[2]][[1]]
  backward_shuffle$Weight <- as.numeric(backward_shuffle$Weight)
  backward_shuffle$edgeIDs <- paste0(backward_shuffle$Node1, '______', backward_shuffle$Sign, '______', backward_shuffle$Node2)
  backward_shuffle <- backward_shuffle[,c(5,4)]
  backward_shuffle <- aggregate(Weight ~ edgeIDs, backward_shuffle, sum)
  backward_shuffle$Weight <- ifelse(backward_shuffle$Weight  > 100, 100, backward_shuffle$Weight)
  
  
  backward_ori <- ori_carni_res[[2]][[1]]
  backward_ori$Weight <- as.numeric(backward_ori$Weight)
  backward_ori$edgeIDs <- paste0(backward_ori$Node1, '______', backward_ori$Sign, '______', backward_ori$Node2)
  backward_ori <- backward_ori[,c(5,4)]
  backward_ori <- aggregate(Weight ~ edgeIDs, backward_ori, sum)
  backward_ori$Weight <- ifelse(backward_ori$Weight  > 100, 100, backward_ori$Weight)
  
  delta_df <- merge(backward_ori, backward_shuffle, by = 'edgeIDs', all = T)
  delta_df[is.na(delta_df)] <- 0
  delta_df$delta <- abs(delta_df$Weight.x - delta_df$Weight.y)
  delta_df$index <- x[[3]]
  
  return(delta_df)
}, ori_carni_res = ori_carni_res)

delta_df_back <- as.data.frame(do.call(rbind, delta_list_backward))
delta_df_back$index <- as.numeric(delta_df_back$index)
delta_df_back <- delta_df_back[order(delta_df_back$index, decreasing = F),]

delta_df_back <- delta_df_back[delta_df_back$index %in% c(2,10,20,30,40,50),]

delta_df_back$index <- factor(delta_df_back$index, levels = unique(delta_df_back$index))

ggplot(delta_df_back, aes(x = index, y = delta, group = index, fill = index)) + 
  geom_boxplot(fill = "lightblue") +
  geom_jitter(alpha = 0.2) +
  stat_summary(fun.y=median, geom="point", shape=23, fill = 'lightgrey', size = 4) +
  theme_minimal()  +
  ggtitle('Backward COSMOS weight comparison') +
  xlab('% of edges shuffled') +
  ylab("Absolute edge weight difference")

ggplot(delta_df_back, aes(x = delta, fill = index)) + 
  geom_density(alpha = 0.33) +
  theme_minimal()

df <- delta_df_for
for(index in unique(df$index))
{
  print(paste0('index : ',index))
  print(median(df[df$index == index,'delta']))
}

df <- delta_df_back
for(index in unique(df$index))
{
  print(paste0('index : ',index))
  print(length(df[df$index == index & df$delta == 0,'delta']) / (length(df[df$index == index & df$delta == 0,'delta']) + length(df[df$index == index & df$delta == 100,'delta'])))
}

df <- delta_df_for
for(index in unique(df$index))
{
  print(paste0('index : ',index))
  print(length(df[df$index == index & df$Weight.y != 0,'delta']))
}

df <- delta_df_back
for(index in unique(df$index))
{
  print(paste0('index : ',index))
  print(length(df[df$index == index & df$Weight.y != 0,'delta']))
}

(149+237+342+340+142+312+195+302+281+288+287+275+277+250)/14

sd(c(149,237,342,340,142,312,195,302,281,288,287,275,277,250))
####full NETWORK

delta_df_for$signed_delta <- delta_df_for$Weight.x - delta_df_for$Weight.y

delta_df_for_X <- delta_df_for[delta_df_for$index == '2',]
length(delta_df_for_X[delta_df_for_X$Weight.x == 0,1])
