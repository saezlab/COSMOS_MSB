library(readr)
library(omicToolsTest)
library(viper)

setwd("~/Dropbox/COSMOS_MSB/")

metabolomic <- as.data.frame(
  read_csv("data/metab_samples_vsn.csv"))

phospho <- as.data.frame(
  read_csv("data/phoshpo_processed.csv"))

omnipath_ptm <- as.data.frame(read_csv("support/omnipath_ptm_20200205.csv"))

omnipath_ptm <- omnipath_ptm[omnipath_ptm$modification %in% c("dephosphorylation","phosphorylation"),]
KSN <- omnipath_ptm[,c(4,3)]
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_type, sep ="_")
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_offset, sep = "")
KSN$sign <- ifelse(omnipath_ptm$modification == "phosphorylation", 1, -1)

KSN_viper <- df_to_viper_regulon(KSN)
KSN_viper <- df_to_viper_regulon(KSN)

row.names(phospho) <- phospho$X1
phospho <- phospho[,-1]
samples <- names(phospho)
phospho <- as.data.frame(t(apply(phospho, 1, scale)))
names(phospho) <- samples


row.names(metabolomic) <- metabolomic$X1
metabolomic <- metabolomic[,-1]

kinases <- as.data.frame(viper(eset = phospho, regulon = KSN_viper, nes = T, minsize = 5, eset.filter = F))
kinases <- kinases[complete.cases(kinases),]
kinases <- kinases[rowSums(kinases) != 0,]

phospho_targets_final <- as.data.frame(
  read_csv("support/phospho_targets_final.csv"))

kinases <- kinases[,names(metabolomic)]
kinases <- kinases[complete.cases(kinases),]
samples <- names(kinases)
metabolomic <- metabolomic[complete.cases(metabolomic),]
metabolomic <- as.data.frame(t(apply(metabolomic, 1, scale)))
names(metabolomic) <- samples

#######
counts_vsn <- as.data.frame(read_csv("data/RNA_counts_vsn.csv"))

dorothea_df <- as.data.frame(
  read_csv("support/DOROTHEA_20200811.csv"))

dorothea_viper <- df_to_viper_regulon(dorothea_df)

row.names(counts_vsn) <- counts_vsn$ID
counts_vsn <- counts_vsn[,-1]
count_sample_names <- names(counts_vsn)

counts_vsn <- as.data.frame(t(apply(counts_vsn, 1, scale)))
names(counts_vsn) <- count_sample_names

TFs <- as.data.frame(viper(eset = counts_vsn, regulon = dorothea_viper, nes = T, minsize = 5, eset.filter = F))
names(TFs) <- gsub("H","KI",names(TFs))
names(TFs) <- gsub("T","TU",names(TFs))

TFs <- TFs[,names(kinases)]
#######
all_features <- as.data.frame(rbind(kinases, TFs))
all_features <- as.data.frame(rbind(all_features, metabolomic))
all_features <- as.data.frame(t(all_features))

corr_mat <- as.data.frame(cor(all_features, method = "spearman"))

edge_list <- list()
k <- 1
for(i in 1:length(corr_mat[,1]))
{
  for( j in i:length(corr_mat[1,]))
  {
    edge_list[[k]] <- c(names(corr_mat)[j], names(corr_mat)[i], corr_mat[i,j])
    k <- k+1
  }
}

edge_df <- as.data.frame(do.call(rbind,edge_list))
names(edge_df) <- c("feature_A","feature_B","correlation")
edge_df$correlation <- as.numeric(as.character(edge_df$correlation))
edge_df_to_write <- edge_df
edge_df <- edge_df[edge_df$correlation != 1,]

top_correlations <- edge_df[abs(edge_df$correlation) > 0.9,]

edge_df_inter <- edge_df[edge_df$feature_A %in% row.names(metabolomic) 
                         & edge_df$feature_B %in% row.names(kinases),]

edge_df_inter <- edge_df_inter[abs(edge_df_inter$correlation) > 0.73,]

reduced_cor_net <- as.data.frame(rbind(top_correlations, edge_df_inter))

write_csv(reduced_cor_net, "results/correlation/reduced_cor_net.csv")
write_csv(edge_df_to_write, "results/correlation/cor_net.csv")

#############################

all_features_tumor <- all_features[grepl("TU",row.names(all_features)),]

all_features_tumor <- all_features_tumor[,rowSums(t(all_features_tumor)) != 0]

corr_mat <- as.data.frame(cor(all_features_tumor, method = "spearman"))

edge_list <- list()
k <- 1
for(i in 1:length(corr_mat[,1]))
{
  for( j in i:length(corr_mat[1,]))
  {
    edge_list[[k]] <- c(names(corr_mat)[j], names(corr_mat)[i], corr_mat[i,j])
    k <- k+1
  }
}

edge_df <- as.data.frame(do.call(rbind,edge_list))
names(edge_df) <- c("feature_A","feature_B","correlation")
edge_df$correlation <- as.numeric(as.character(edge_df$correlation))
edge_df_to_write <- edge_df
edge_df <- edge_df[edge_df$correlation != 1,]

top_correlations <- edge_df[abs(edge_df$correlation) > 0.97,]

edge_df_inter <- edge_df[edge_df$feature_A %in% row.names(metabolomic) 
                         & edge_df$feature_B %in% row.names(kinases),]

edge_df_inter <- edge_df_inter[abs(edge_df_inter$correlation) > 0.92,]

reduced_cor_net <- as.data.frame(rbind(top_correlations, edge_df_inter))

write_csv(reduced_cor_net, "results/correlation/revision_reduced_cor_net_tumor.csv")
write_csv(edge_df_to_write, "results/correlation/revision_cor_net_tumor.csv")
# 
# ###############################
# 
# all_features_healthy <- all_features[grepl("KI",row.names(all_features)),]
# 
# corr_mat <- as.data.frame(cor(all_features_healthy, method = "spearman"))
# 
# edge_list <- list()
# k <- 1
# for(i in 1:length(corr_mat[,1]))
# {
#   for( j in i:length(corr_mat[1,]))
#   {
#     edge_list[[k]] <- c(names(corr_mat)[j], names(corr_mat)[i], corr_mat[i,j])
#     k <- k+1
#   }
# }
# 
# edge_df <- as.data.frame(do.call(rbind,edge_list))
# names(edge_df) <- c("feature_A","feature_B","correlation")
# edge_df$correlation <- as.numeric(as.character(edge_df$correlation))
# edge_df_to_write <- edge_df
# edge_df <- edge_df[edge_df$correlation != 1,]
# 
# top_correlations <- edge_df[abs(edge_df$correlation) > 0.97,]
# 
# edge_df_inter <- edge_df[edge_df$feature_A %in% row.names(metabolomic) 
#                          & edge_df$feature_B %in% row.names(kinases),]
# 
# edge_df_inter <- edge_df_inter[abs(edge_df_inter$correlation) > 0.92,]
# 
# reduced_cor_net <- as.data.frame(rbind(top_correlations, edge_df_inter))
# 
# write_csv(reduced_cor_net, "~/Dropbox/kidney_cancer_multiomic_pipeline/results/multi_omic/final_sample/correlation/reduced_cor_net_healthy.csv")
# write_csv(edge_df_to_write, "~/Dropbox/kidney_cancer_multiomic_pipeline/results/multi_omic/final_sample/correlation/cor_net_healthy.csv")
