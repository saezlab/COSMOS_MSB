library(readr)
library(limma)
library(viper)
library('org.Hs.eg.db')

setwd("~/Dropbox/COSMOS_MSB//")
source("scripts/limma_functions.R")
source("scripts/viper_functions.R")

phospho <- as.data.frame(
  read_delim("data/CPTAC_phosphosites.txt", "\t", escape_double = FALSE, trim_ws = TRUE))
row.names(phospho) <- phospho$psite
phospho <- phospho[,-1]


targets <- as.data.frame(
  read_delim("support//CPTAC_samples.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE))

names(targets) <- c("sample","condition")

unique(targets$condition)

limmaRes <- runLimma(phospho, targets, comparisons = list(c(1,-2)))

ttop <- ttopFormatter(topTable(limmaRes[[1]], coef = 1, number = 20284, adjust.method = "fdr"))

write_csv(ttop,"data/CPTAC_phospho_ttop.csv")

omnipath_ptm <- as.data.frame(read_csv("support/omnipath_ptm_20200205.csv"))

omnipath_ptm <- omnipath_ptm[omnipath_ptm$modification %in% c("dephosphorylation","phosphorylation"),]
KSN <- omnipath_ptm[,c(4,3)]
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_type, sep ="_")
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_offset, sep = "")
KSN$sign <- ifelse(omnipath_ptm$modification == "phosphorylation", 1, -1)

KSN_viper <- df_to_viper_regulon(KSN)

#Prepare the measurments for viper
viper_expression <- ttop$t
names(viper_expression) <- ttop$ID

Kin_activity <- as.data.frame(viper(eset = viper_expression, regulon = KSN_viper, minsize = 5, adaptive.size = F, eset.filter = F))
kinases <- row.names(Kin_activity) 

write.csv(Kin_activity,"data/CPTAC_Kinase_Activities_omnipath.csv")

top_kinase <- Kin_activity[abs(Kin_activity$V1) > 1.7,,drop = F]
top_kinase <- as.data.frame(t(top_kinase))


mapping_symbole_to_entrez <- mapIds(org.Hs.eg.db, names(top_kinase), 'ENTREZID', 'SYMBOL')
for(i in 1:length(names(top_kinase)))
{
  names(top_kinase)[i] <- mapping_symbole_to_entrez[names(top_kinase)[i]]
}
names(top_kinase) <- paste0("X",names(top_kinase))

write_csv(top_kinase,"data/CPTAC_top_kinase_carni_input.csv")
