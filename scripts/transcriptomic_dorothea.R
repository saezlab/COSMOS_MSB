library(dorothea)
library(readr)
library(viper)

# dorothea_df<- as.data.frame(dorothea_hs[dorothea_hs$confidence %in% c("A","B","C","D"),c(3,1,4)])

setwd("~/Dropbox/COSMOS_MSB/")

# write_csv(dorothea_df, "support/DOROTHEA_20200811.csv")
dorothea_df <- as.data.frame(
  read_csv("support/DOROTHEA_20200811.csv"))

source("scripts/limma_functions.R")
source("scripts/viper_functions.R")

dorothea_viper <- df_to_viper_regulon(dorothea_df)

ttop_tumorvshealthy <- as.data.frame(
  read_csv("data/RNA_ttop_tumorvshealthy.csv"))

RNAseq_entrez_to_symbol <- as.data.frame(read_delim("support/RNAseq_entrez_to_symbol", 
                                                    "\t", escape_double = FALSE, col_types = cols(`yourlist:M20191127216DA2B77BFBD2E6699CA9B6D1C41EB259129CL` = col_character()), 
                                                    trim_ws = TRUE)) #from uniprot 20191127

names(RNAseq_entrez_to_symbol)[1] <- "ID"

ttop_tumorvshealthy <- merge(ttop_tumorvshealthy, RNAseq_entrez_to_symbol[,c(1,6)])
ttop_tumorvshealthy <- ttop_tumorvshealthy[,c(8,2:7)]
names(ttop_tumorvshealthy)[1] <- "ID"
ttop_tumorvshealthy$ID <- gsub(" .*","",ttop_tumorvshealthy$ID)
ttop_tumorvshealthy <- unique(ttop_tumorvshealthy)

eset <- ttop_tumorvshealthy$t
names(eset) <- ttop_tumorvshealthy$ID

viperRes <- as.data.frame(viper(eset = eset, regulon = dorothea_viper, minsize = 10, adaptive.size = F, eset.filter = F, pleiotropy = T))
viperRes$TF <- row.names(viperRes)

viperRes <- viperRes[,c(2,1)]
names(viperRes) <- c("ID","NES")

# TF_scores_prerevision <- as.data.frame(read_csv("results/transcriptomic/TF_scores.csv"))

write_csv(viperRes, "data/RNA_TF_scores.csv")

# comp <- merge(viperRes, TF_scores_prerevision, by = "ID", all = T)
