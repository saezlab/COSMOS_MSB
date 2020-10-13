library(readr)
library(vsn)
library(limma)
library(viper)
library(dorothea)
library('org.Hs.eg.db')

setwd("~/Dropbox/COSMOS_MSB//")
source("scripts/limma_functions.R")
source("scripts/viper_functions.R")

RNA_rpkm_tumor_normal <- as.data.frame(
  read_delim("data/CPTAC_RNA_rpkm_tumor_normal.tsv", 
             "\t", escape_double = FALSE, trim_ws = TRUE))

doublons <- RNA_rpkm_tumor_normal[duplicated(RNA_rpkm_tumor_normal$geneID),1]
RNA_rpkm_tumor_normal <- RNA_rpkm_tumor_normal[!(RNA_rpkm_tumor_normal$geneID %in% doublons),]

row.names(RNA_rpkm_tumor_normal) <- RNA_rpkm_tumor_normal$geneID

RNA_rpkm_tumor_normal <- RNA_rpkm_tumor_normal[,-1]

RNA_clinical <- as.data.frame(
  read_csv("support//CPTAC_RNA_clinical.csv", 
                         col_names = FALSE))

sample_mapping <- RNA_clinical$X1
names(sample_mapping) <- RNA_clinical$X3

RNA_rpkm_tumor_normal <- RNA_rpkm_tumor_normal[,names(RNA_rpkm_tumor_normal) %in% RNA_clinical$X3]

names(RNA_rpkm_tumor_normal) <- sapply(names(RNA_rpkm_tumor_normal), function(x, sample_mapping)
  {
  return(sample_mapping[x])
},sample_mapping = sample_mapping, simplify = T)

plot(hist(as.matrix(log2(RNA_rpkm_tumor_normal)),breaks = 10000))

RNA_rpkm_tumor_normal[log2(RNA_rpkm_tumor_normal) <= 0] <- NA

RNA_rpkm_tumor_normal <- RNA_rpkm_tumor_normal[rowSums(is.na(RNA_rpkm_tumor_normal)) < 170,]

fit <- vsnMatrix(as.matrix(RNA_rpkm_tumor_normal))
meanSdPlot(fit)
RNA_rpkm_tumor_normal_vsn <- as.data.frame(vsn::predict(fit,as.matrix(RNA_rpkm_tumor_normal)))

targets <- RNA_clinical[,c(1,2)]
names(targets) <- c("sample","condition")

unique(targets$condition)

limmaRes <- runLimma(RNA_rpkm_tumor_normal_vsn, targets, comparisons = list(c(2,-1)))

ttop <- ttopFormatter(topTable(limmaRes[[1]], coef = 1, number = 14935, adjust.method = "fdr"))

########## DOROTHEA

dorothea_df <- as.data.frame(
  read_csv("support/DOROTHEA_20200811.csv"))

dorothea_viper <- df_to_viper_regulon(dorothea_df)

eset <- ttop$t
names(eset) <- ttop$ID

viperRes <- as.data.frame(viper(eset = eset, regulon = dorothea_viper, minsize = 10, adaptive.size = F, eset.filter = F, pleiotropy = T))
viperRes$TF <- row.names(viperRes)

viperRes <- viperRes[,c(2,1)]
names(viperRes) <- c("ID","NES")

viperRes <- viperRes[abs(viperRes$NES) > 1.7,]

TF_inputs <- as.data.frame(t(viperRes[,2,drop = F]))

mapping_symbole_to_entrez <- mapIds(org.Hs.eg.db, names(TF_inputs), 'ENTREZID', 'SYMBOL')
for(i in 1:length(names(TF_inputs)))
{
  names(TF_inputs)[i] <- mapping_symbole_to_entrez[names(TF_inputs)[i]]
}

names(TF_inputs) <- paste("X", names(TF_inputs), sep = "")

write_csv(TF_inputs, "data/CPTAC_carni_TF_inputs.csv")

mapping_symbole_to_entrez <- mapIds(org.Hs.eg.db, ttop$ID, 'ENTREZID', 'SYMBOL')
for(i in 1:length(ttop$ID))
{
  ttop[i,"ID"] <- mapping_symbole_to_entrez[ttop[i,"ID"]]
}
ttop$ID <- paste0("X",ttop$ID)

write_csv(ttop,"data/CPTAC_RNA_ttop_symbol.csv")
