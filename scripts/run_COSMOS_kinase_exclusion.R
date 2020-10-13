library(readr)
library(igraph)
library(CARNIVAL)
library(biomaRt)

setwd("~/Dropbox/COSMOS_MSB/")
source("scripts/COSMOS_functions.R")

my_stat <- "t" 
my_threshold <- 1

##Dorothea/viper
load("support/dorothea_pkg_df.RData")

#Load inputs
meta_network <- as.data.frame(
  read_csv("support/metaPKN_filtered.csv"))

metab_input_carnival <- as.data.frame(read_csv("data/metab_input_COSMOS.csv"))

signaling_input_carnival <- as.data.frame(read_csv("data/signaling_input_COSMOS.csv"))

print(names(metab_input_carnival))
print(names(signaling_input_carnival))

ttop_RNA <- as.data.frame(read_csv("data/RNA_ttop_tumorvshealthy.csv"))
ttop_RNA$ID <- paste0("X",ttop_RNA$ID)

#Preprocess pkn

signaling_input_carnival <- filter_inputs(signaling_input_carnival, meta_network)

#Exclude Kinases 

signaling_input_carnival <- signaling_input_carnival[,(names(signaling_input_carnival) %in% dorothea_df$source_genesymbol)]

metab_input_carnival <- filter_inputs(metab_input_carnival, meta_network)

meta_network <- downstream_neighbours(meta_network, 8, names(signaling_input_carnival))
###############

meta_network <- filter_TF_sign(meta_network = meta_network, ttop_RNA = ttop_RNA, inputs = signaling_input_carnival, TF_targets = dorothea_df)

setwd("~/Dropbox/COSMOS_MSB/results/COSMOS_result_kinase_exclusion/")
###############

CARNIVAL_Result <- runCARNIVAL(inputObj = sign(signaling_input_carnival),
                               measObj = metab_input_carnival,
                               netObj = meta_network,
                               solverPath = "~/Documents/cplex",
                               solver = "cplex",
                               timelimit = 3600,
                               mipGAP = 0.2)

save.image("carni_forward_prerun_kinase_exclusion.RData")

TF_signs <- get_TF_sign_from_CARNI(CARNIVAL_Result, dorothea_df)

meta_network <- filter_TF_sign(meta_network = meta_network, ttop_RNA = ttop_RNA, inputs = TF_signs, TF_targets = dorothea_df)

CARNIVAL_Result_rerun <- runCARNIVAL(inputObj = sign(signaling_input_carnival),
                                     measObj = metab_input_carnival,
                                     netObj = meta_network,
                                     solverPath = "~/Documents/cplex",
                                     solver = "cplex",
                                     timelimit = 7200,
                                     mipGAP = 0.2)

save.image("carni_forwardrun_kinase_exclusion.RData")

meta_network <- as.data.frame(
  read_csv("../../support/metaPKN_filtered.csv"))


meta_network <- downstream_neighbours(meta_network, 8, names(metab_input_carnival))

######
######
######
meta_network <- filter_TF_sign(meta_network = meta_network, ttop_RNA = ttop_RNA, inputs = signaling_input_carnival, TF_targets = dorothea_df)

meta_network <- filter_TF_sign(meta_network = meta_network, ttop_RNA = ttop_RNA, inputs = TF_signs, TF_targets = dorothea_df)

####
####
####

CARNIVAL_Result_2 <- runCARNIVAL(inputObj = sign(metab_input_carnival),
                                 measObj = signaling_input_carnival,
                                 netObj = meta_network,
                                 solverPath = "~/Documents/cplex",
                                 solver = "cplex",
                                 timelimit = 10800,
                                 mipGAP = 0.2)

save.image("carni_doublerun_kinase_exclusion.RData")