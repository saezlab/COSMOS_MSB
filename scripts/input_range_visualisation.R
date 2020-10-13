library(pheatmap)

setwd("~/Dropbox/COSMOS_MSB/")

load("data/input_range_carnival_input_list.Rdata")

sig_1 <- as.data.frame(t(carnival_input_list[[1]][[1]]))
sig_2 <- as.data.frame(t(carnival_input_list[[2]][[1]]))

sig_1$ID <- row.names(sig_1)
sig_2$ID <- row.names(sig_2)

sig_hm <- merge(sig_1, sig_2, by = "ID", all = T)

for(i in 3:25)
{
  sig_i <- as.data.frame(t(carnival_input_list[[i]][[1]]))
  sig_i$ID <- row.names(sig_i)
  
  sig_hm <- merge(sig_hm, sig_i, by = "ID", all = T)
}

row.names(sig_hm) <- sig_hm$ID
sig_hm <- sig_hm[,-1]

names(sig_hm) <- paste0("threshold_",1:25)
sig_hm[!is.na(sig_hm)] <- 1
sig_hm[is.na(sig_hm)] <- 0

order_rows <- rowSums(sig_hm)
order_rows <- order(order_rows, decreasing = F)

sig_hm <- sig_hm[order_rows,]

pheatmap(sig_hm, cluster_cols = F, cluster_rows = F, show_rownames = F, color = c("white","royalblue"), legend = F, border_color = "white")

met_1 <- as.data.frame(t(carnival_input_list[[1]][[2]]))
met_2 <- as.data.frame(t(carnival_input_list[[2]][[2]]))

met_1$ID <- row.names(met_1)
met_2$ID <- row.names(met_2)

met_hm <- merge(met_1, met_2, by = "ID", all = T)

for(i in 3:25)
{
  met_i <- as.data.frame(t(carnival_input_list[[i]][[2]]))
  met_i$ID <- row.names(met_i)
  
  met_hm <- merge(met_hm, met_i, by = "ID", all = T)
}

met_hm$ID <- gsub("___[a-z]____","",met_hm$ID)
met_hm <- met_hm[!duplicated(met_hm$ID),]

row.names(met_hm) <- met_hm$ID
met_hm <- met_hm[,-1]

names(met_hm) <- paste0("threshold_",1:25)
met_hm[!is.na(met_hm)] <- 1
met_hm[is.na(met_hm)] <- 0


order_rows <- rowSums(met_hm)
order_rows <- order(order_rows, decreasing = F)

met_hm <- met_hm[order_rows,]

pheatmap(met_hm, cluster_cols = F, cluster_rows = F, show_rownames = F, color = c("white","royalblue"), legend = F, border_color = "white")
