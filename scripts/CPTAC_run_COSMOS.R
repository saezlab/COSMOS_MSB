library(readr)
library(dplyr)
library(reshape2)
library(biomaRt)
library('org.Hs.eg.db')
library(igraph)
library(CARNIVAL)
source("/home/ad234505/COSMOS_revisions/CCRCC_CPTAC/revision_COSMOS_functions.R")

top_kinases <- as.data.frame(read_csv("/home/ad234505/COSMOS_revisions/CCRCC_CPTAC/ccrcc_CPTAC_top_kinase_carni_input.csv"))

top_TF <- as.data.frame(read_csv("/home/ad234505/COSMOS_revisions/CCRCC_CPTAC/carni_TF_inputs.csv"))

ttop_RNA <- as.data.frame(
  read_csv("/home/ad234505/COSMOS_revisions/CCRCC_CPTAC/ttop_CPTAC_CCRCC.csv")
)

meta_network <- as.data.frame(
  read_csv("/home/ad234505/COSMOS_revisions/CCRCC_CPTAC/meta_network_CCRCC_CPTAC_expfiltered.csv"))

load("/home/ad234505/COSMOS_revisions/CCRCC_CPTAC/dorothea_pkg_df.RData") #ABCD

########## COSMOS #################


setwd("/home/ad234505/COSMOS_revisions/CCRCC_CPTAC/")

my_stat <- "t" 
my_threshold <- 1

#Preprocess pkn

top_kinases <- filter_inputs(top_kinases, meta_network)

top_TF <- filter_inputs(top_TF, meta_network)

meta_network <- downstream_neighbours(meta_network, 8, names(top_kinases))
###############

meta_network <- filter_TF_sign(meta_network = meta_network, ttop_RNA = ttop_RNA, inputs = top_kinases, TF_targets = dorothea_df)
###############

CARNIVAL_Result <- runCARNIVAL(inputObj = sign(top_kinases),
                               measObj = top_TF,
                               netObj = meta_network,
                               solverPath = "cplex",
                               solver = "cplex",
                               timelimit = 3600,
                               mipGAP = 0.2)

save.image("carni_forwardrun_res_fullomni_fullfilter_first_newDoro_long.RData")

TF_signs <- get_TF_sign_from_CARNI(CARNIVAL_Result, dorothea_df)

meta_network <- filter_TF_sign(meta_network = meta_network, ttop_RNA = ttop_RNA, inputs = TF_signs, TF_targets = dorothea_df)

CARNIVAL_Result_rerun <- runCARNIVAL(inputObj = sign(top_kinases),
                                     measObj = top_TF,
                                     netObj = meta_network,
                                     solverPath = "cplex",
                                     solver = "cplex",
                                     timelimit = 7200,
                                     mipGAP = 0.2)

save.image("carni_forwardrun_res_fullomni_fullfilter_second_newDoro_long.RData")

meta_network <- as.data.frame(
  read_csv("/home/ad234505/COSMOS_revisions/CCRCC_CPTAC/meta_network_CCRCC_CPTAC_expfiltered.csv"))


meta_network <- downstream_neighbours(meta_network, 8, names(top_TF))

######
######
######
meta_network <- filter_TF_sign(meta_network = meta_network, ttop_RNA = ttop_RNA, inputs = top_kinases, TF_targets = dorothea_df)

meta_network <- filter_TF_sign(meta_network = meta_network, ttop_RNA = ttop_RNA, inputs = TF_signs, TF_targets = dorothea_df)

####
####
####

CARNIVAL_Result_2 <- runCARNIVAL(inputObj = sign(top_TF),
                                 measObj = top_kinases,
                                 netObj = meta_network,
                                 solverPath = "cplex",
                                 solver = "cplex",
                                 timelimit = 21600,
                                 mipGAP = 0.2)

save.image("carni_doublerun_res_fullomni_met_expfiltered_newDoro_long.RData")