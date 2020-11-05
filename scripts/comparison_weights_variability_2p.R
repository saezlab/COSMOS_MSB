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
library(ggplot2)

setwd("~/Dropbox/COSMOS_MSB/results/variability/")
source("../../scripts/format_cosmos_res_function.R")
source("../../scripts/COSMOS_functions.R")
source("../../scripts/piano_function.R")

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

files <- list.files("../shuffler/2/", recursive = T)
files <- files[grepl("doublerun_res",files)]
files <- paste0("../shuffler/2/",files)

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
  run_index <- gsub("[/].*","shuffle",file)
  run_index <- gsub("[.][.]","2",run_index)
  
  carni_res_list[[i]] <- list(cosmos_forward, cosmos_backward)
  names(carni_res_list)[i] <- run_index
  i <- i+1
}

for(i in 1:(length(carni_res_list)-1))
{
  if(i == 1)
  {
    sif_combined <-  carni_res_list[[i]][[1]][[1]]
    sif_combined$edgeIDs <- paste0(sif_combined$Node1, '______', sif_combined$Sign, '______', sif_combined$Node2)
    sif_combined <- sif_combined[,c(5,4)]
    sif_combined$Weight <- as.numeric(sif_combined$Weight)
    sif_combined <- aggregate(Weight ~ edgeIDs, sif_combined, sum)
  }
  else
  {
    print(i)
    new_sif <- carni_res_list[[i]][[1]][[1]]
    new_sif$edgeIDs <- paste0(new_sif$Node1, '______', new_sif$Sign, '______', new_sif$Node2)
    new_sif <- new_sif[,c(5,4)]
    new_sif$Weight <- as.numeric(new_sif$Weight)
    new_sif <- aggregate(Weight ~ edgeIDs, new_sif, sum)
    
    sif_combined <- merge(sif_combined, new_sif, by = 'edgeIDs', all = T)
    sif_combined[is.na(sif_combined)] <- 0
    names(sif_combined) <- c('edgeIDs', paste0('run_',1:i))
  }
}

sif_combined$Weight <- rowMeans(sif_combined[,-1])
sif_combined <- sif_combined[,c(1,length(sif_combined[1,]))]

full <- carni_res_list[[1]][[1]][[1]]
p2 <- carni_res_list[[25]][[1]][[1]]

full$edgeIDs <- paste0(full$Node1, '______', full$Sign, '______', full$Node2)
p2$edgeIDs <- paste0(p2$Node1, '______', p2$Sign, '______', p2$Node2)
p2 <- p2[,c(5,4)]
p2$Weight <- as.numeric(p2$Weight)
p2 <- aggregate(Weight ~ edgeIDs, p2, sum)

comparison <- merge(full[,c(5,4)], p2, by = 'edgeIDs', all = T)
comparison[is.na(comparison)] <- 0

comparison$delta  <- abs(as.numeric(comparison$Weight.x) - as.numeric(comparison$Weight.y))
median_delta <- median(comparison$delta)

ggplot(comparison, aes(x = delta)) +
  geom_density(fill = 'lightblue') +
  geom_vline(xintercept = median_delta) +
  theme_minimal()

comparison_combined <- merge(sif_combined, p2, by = 'edgeIDs', all = T)
comparison_combined[is.na(comparison_combined)] <- 0
comparison_combined$delta  <- abs(as.numeric(comparison_combined$Weight.x) - as.numeric(comparison_combined$Weight.y))
median_delta <- median(comparison_combined$delta, na.rm = T)

ggplot(comparison_combined, aes(x = delta)) +
  geom_density(fill = 'lightblue') +
  geom_vline(xintercept = median_delta) +
  theme_minimal()
