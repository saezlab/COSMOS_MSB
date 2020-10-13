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
library(visNetwork)

setwd("~/Dropbox/COSMOS_MSB/")

CARNIVAL_Result_rerun <- readRDS("results/forward_run.rds")
CARNIVAL_Result_2 <- readRDS("results/backward_run.rds")

metab_input <- as.data.frame(
  read_csv("data/metab_input_COSMOS.csv"))

signaling_input <- as.data.frame(read_csv("data/signaling_input_COSMOS.csv"))

source("scripts/format_cosmos_res_function.R")

metab_to_pubchem <- as.data.frame(read_csv("support/metab_to_pubchem.csv"))

cosmos_forward <- format_COSMOS_res(cosmos_res = CARNIVAL_Result_rerun,
                                    metab_mapping = metab_to_pubchem,
                                    measured_nodes = c(names(metab_input),names(signaling_input)),
                                    omnipath_ptm = read_rds("support/omnipath_ptm.rds"), gene_mapping = "else")

cosmos_backward <- format_COSMOS_res(cosmos_res = CARNIVAL_Result_2,
                                    metab_mapping = metab_to_pubchem,
                                    measured_nodes = c(names(metab_input),names(signaling_input)),
                                    omnipath_ptm = read_rds("support/omnipath_ptm.rds"), gene_mapping = "else")

write_csv(cosmos_forward[[1]], "sif_first_newDoro_long.csv")

write_csv(cosmos_forward[[2]], "att_first_newDoro_long.csv")

write_csv(cosmos_backward[[1]], "sif_sec_newDoro_long.csv")

write_csv(cosmos_backward[[2]], "att_sec_newDoro_long.csv")

full_sif <- as.data.frame(unique(rbind(cosmos_forward[[1]],cosmos_backward[[1]])))
full_att <- as.data.frame(unique(rbind(cosmos_forward[[2]][,-6],cosmos_backward[[2]][,-6])))
full_att <- full_att[full_att$AvgAct != 0,]
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

write_csv(full_sif, "subnet_combined_sif_full_newDoro_long.csv")

write_csv(full_att, "subnet_combined_att_full_newDoro_long.csv")
