library(readr)
library(dplyr)
library(reshape2)
library(biomaRt)
library('org.Hs.eg.db')
library(igraph)
library(CARNIVAL)
source("/home/ad234505/COSMOS_revisions/CCRCC_CPTAC/revision_COSMOS_functions.R")


kinases_CCRCC_CPTAC <- as.data.frame(
  read_delim("/home/ad234505/COSMOS_revisions/CCRCC_CPTAC/CCRCC_CPTAC_Kinase_Activities.txt", 
             "\t", escape_double = FALSE, trim_ws = TRUE))

kinases_CCRCC_CPTAC <- kinases_CCRCC_CPTAC[,c(2,3,5)]
kinases_CCRCC_CPTAC <- dcast(kinases_CCRCC_CPTAC, formula = kinase~sample)
row.names(kinases_CCRCC_CPTAC) <- kinases_CCRCC_CPTAC$kinase
kinases_CCRCC_CPTAC <- kinases_CCRCC_CPTAC[,-1]
kinases_CCRCC_CPTAC <- kinases_CCRCC_CPTAC[rowSums(is.na(kinases_CCRCC_CPTAC)) <= length(kinases_CCRCC_CPTAC[1,])*90/100,]

kinase_meanoversd <- apply(kinases_CCRCC_CPTAC,1,function(x){mosd <- mean(x,na.rm = T) / sd(x, na.rm = T)})


plot(density(kinase_meanoversd))

bot_10pcent_kinase <- sort(kinase_meanoversd)[1:(length(kinase_meanoversd)*0.1)]
top_10pcent_kinase <- sort(kinase_meanoversd, decreasing = T)[1:(length(kinase_meanoversd)*0.1)] 

top_kinases <- c(bot_10pcent_kinase, top_10pcent_kinase)


top_kinases <- as.data.frame(t(top_kinases))

symbols <- names(top_kinases)

# use mapIds method to obtain Entrez IDs
mapping_symbole_to_entrez <- mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
for(i in 1:length(names(top_kinases)))
{
  names(top_kinases)[i] <- mapping_symbole_to_entrez[names(top_kinases)[i]]
}


names(top_kinases) <- paste("X", names(top_kinases), sep = "")

top_TF <- as.data.frame(read_csv("/home/ad234505/COSMOS_revisions/CCRCC_CPTAC/carni_TF_inputs.csv"))

ttop_RNA <- as.data.frame(
  read_csv("/home/ad234505/COSMOS_revisions/CCRCC_CPTAC/ttop_CPTAC_CCRCC.csv")
)

meta_network <- as.data.frame(
  read_csv("/home/ad234505/COSMOS_revisions/CCRCC_CPTAC/meta_network_CCRCC_CPTAC_expfiltered.csv"))

load("/home/ad234505/COSMOS_revisions/CCRCC_CPTAC/dorothea_pkg_df.RData") #ABCD

########## END OF PREPROCESSING ##################

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