library(pheatmap)
library(readr)

setwd("~/Dropbox/COSMOS_MSB/")

signaling_input_COSMOS <- as.data.frame(read_csv("data/signaling_input_COSMOS.csv"))
signaling_input_COSMOS <- as.data.frame(t(signaling_input_COSMOS))
signaling_input_COSMOS$ID <- row.names(signaling_input_COSMOS)
names(signaling_input_COSMOS)[1] <- 'original'

metab_input_COSMOS <- as.data.frame(read_csv("data/metab_input_COSMOS.csv"))
metab_input_COSMOS <- as.data.frame(t(metab_input_COSMOS))
metab_input_COSMOS$ID <- row.names(metab_input_COSMOS)
names(metab_input_COSMOS)[1] <- 'original'
metab_input_COSMOS <- metab_input_COSMOS[row.names(metab_input_COSMOS) %in% row.names(as.data.frame(t(carnival_input_list[[1]][[2]])))]

load("data/input_range_carnival_input_list.Rdata")
metab_input_COSMOS <- metab_input_COSMOS[row.names(metab_input_COSMOS) %in% row.names(as.data.frame(t(carnival_input_list[[1]][[2]]))),]

sig_1 <- as.data.frame(t(carnival_input_list[[1]][[1]]))
sig_2 <- as.data.frame(t(carnival_input_list[[2]][[1]]))

sig_1$ID <- row.names(sig_1)
sig_2$ID <- row.names(sig_2)

sig_hm <- merge(sig_1, sig_2, by = "ID", all = T)

ori_included <- F
for(i in 3:25)
{
  sig_i <- as.data.frame(t(carnival_input_list[[i]][[1]]))
  sig_i$ID <- row.names(sig_i)
  
  if(length(sig_i[,1]) > length(signaling_input_COSMOS[,1])  | ori_included)
  {
    sig_hm <- merge(sig_hm, sig_i, by = "ID", all = T)
  } else
  {
    sig_hm <- merge(sig_hm, signaling_input_COSMOS, by = "ID", all = T)
    ori_included <- T
    sig_hm <- merge(sig_hm, sig_i, by = "ID", all = T)
  }
}

row.names(sig_hm) <- sig_hm$ID
sig_hm <- sig_hm[,-1]

names(sig_hm) <- gsub('^V.*','threshold',names(sig_hm))
j <- 1
for(i in 1:length(names(sig_hm)))
{
  names(sig_hm)[i] <- gsub('threshold',paste0('threshold_',j), names(sig_hm)[i])
  if(names(sig_hm)[i] != 'original')
  {
    j <- j+1
  }
}
# names(sig_hm) <- paste0(names(sig_hm),'_',1:26)
sig_hm[!is.na(sig_hm)] <- 1
sig_hm[is.na(sig_hm)] <- 0
sig_hm$original <- ifelse(sig_hm$original == 1, 2, 0)

order_rows <- rowSums(sig_hm)
order_rows <- order(order_rows, decreasing = F)

sig_hm <- sig_hm[order_rows,]

pheatmap(sig_hm, cluster_cols = F, cluster_rows = F, show_rownames = F, color = c("white","royalblue",'red'), legend = F, border_color = "white")

met_1 <- as.data.frame(t(carnival_input_list[[1]][[2]]))
met_2 <- as.data.frame(t(carnival_input_list[[2]][[2]]))

met_1$ID <- row.names(met_1)
met_2$ID <- row.names(met_2)

met_hm <- merge(met_1, met_2, by = "ID", all = T)

ori_included <- F
for(i in 3:25)
{
  met_i <- as.data.frame(t(carnival_input_list[[i]][[2]]))
  met_i$ID <- row.names(met_i)
  
  if(length(met_i[,1]) > length(metab_input_COSMOS[,1])  | ori_included)
  {
    met_hm <- merge(met_hm, met_i, by = "ID", all = T)
  } else
  {
    met_hm <- merge(met_hm, metab_input_COSMOS, by = "ID", all = T)
    ori_included <- T
    met_hm <- merge(met_hm, met_i, by = "ID", all = T)
  }
}

row.names(met_hm) <- met_hm$ID
met_hm <- met_hm[,-1]

names(met_hm) <- gsub('^V.*','threshold',names(met_hm))
j <- 1
for(i in 1:length(names(met_hm)))
{
  names(met_hm)[i] <- gsub('threshold',paste0('threshold_',j), names(met_hm)[i])
  if(names(met_hm)[i] != 'original')
  {
    j <- j+1
  }
}
# names(met_hm) <- paste0(names(met_hm),'_',1:26)
met_hm[!is.na(met_hm)] <- 1
met_hm[is.na(met_hm)] <- 0
met_hm$original <- ifelse(met_hm$original == 1, 2, 0)

order_rows <- rowSums(met_hm)
order_rows <- order(order_rows, decreasing = F)

met_hm <- met_hm[order_rows,]

pheatmap(met_hm, cluster_cols = F, cluster_rows = F, show_rownames = F, color = c("white","royalblue",'red'), legend = F, border_color = "white")
