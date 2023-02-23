library(ggplot2)

setwd("~/Dropbox/COSMOS_MSB/results/shuffler/")

load("../COSMOS_result/COSMOS_res_session.RData")

original_res_for <- as.data.frame(CARNIVAL_Result_rerun$weightedSIF)
original_res_back <- as.data.frame(CARNIVAL_Result_2$weightedSIF)
original_PKN <- meta_network
original_PKN$edgeIDs <- paste0(original_PKN$source, '______', original_PKN$interaction, '______', original_PKN$target)

files <- list.files(".", recursive = T)
files <- files[grepl("doublerun_res",files)]

delta_forward_df_list <- lapply(files, function(x, original_PKN, original_res_for){
  print(x)
  load(x)
  
  shuffled_res_for <- as.data.frame(CARNIVAL_Result_rerun$weightedSIF)
  shuffled_PKN <- meta_network
  shuffled_PKN$edgeIDs <- paste0(shuffled_PKN$source, '______', shuffled_PKN$interaction, '______', shuffled_PKN$target)
  
  # common_PKN <- merge(original_PKN, shuffled_PKN, by = 'edgeIDs', all.x = T)
  # common_PKN <- common_PKN[,1,drop = F]
  common_PKN <- original_PKN[,4,drop = F]
  
  original_res_for$edgeIDs <- paste0(original_res_for$Node1, '______', original_res_for$Sign, '______', original_res_for$Node2)
  original_res_for <- original_res_for[,c(5,4)]
  
  shuffled_res_for$edgeIDs <- paste0(shuffled_res_for$Node1, '______', shuffled_res_for$Sign, '______', shuffled_res_for$Node2)
  shuffled_res_for <- shuffled_res_for[,c(5,4)]
  
  common_PKN <- merge(common_PKN, original_res_for, by = 'edgeIDs', all.x = T)
  common_PKN$Weight <- ifelse(is.na(common_PKN$Weight), 0, common_PKN$Weight)
  
  common_PKN <- merge(common_PKN, shuffled_res_for, by = 'edgeIDs', all.x = T)
  common_PKN$Weight.y <- ifelse(is.na(common_PKN$Weight.y), 0, common_PKN$Weight.y)
  
  common_PKN$Weight.x <- as.numeric(common_PKN$Weight.x)
  common_PKN$Weight.y <- as.numeric(common_PKN$Weight.y)
  
  common_PKN$delta <- abs(common_PKN$Weight.x - common_PKN$Weight.y)
  common_PKN$percent <- gsub('/.*','',x)
  
  common_PKN$delta <- ifelse(common_PKN$edgeIDs %in% shuffled_PKN$edgeIDs, common_PKN$delta, 100)
  
  return(common_PKN)
}, original_PKN = original_PKN, original_res_for = original_res_for)

full_forward_delta_df <- as.data.frame(do.call(rbind, delta_forward_df_list))

sub_full_forward_delta_df <- full_forward_delta_df[full_forward_delta_df$percent %in% c(2,10,20,30,40,50),]

ggplot(sub_full_forward_delta_df, aes(x = delta, fill = percent)) + 
  geom_density(alpha = 0.33) +
  theme_minimal()

delta_backward_df_list <- lapply(files, function(x, original_PKN, original_res_back){
  print(x)
  load(x)
  
  shuffled_res_back <- as.data.frame(CARNIVAL_Result_2$weightedSIF)
  shuffled_PKN <- meta_network
  shuffled_PKN$edgeIDs <- paste0(shuffled_PKN$source, '______', shuffled_PKN$interaction, '______', shuffled_PKN$target)
  
  # common_PKN <- merge(original_PKN, shuffled_PKN, by = 'edgeIDs', all.x = T)
  # common_PKN <- common_PKN[,1,drop = F]
  common_PKN <- original_PKN[,4,drop = F]
  
  original_res_back$edgeIDs <- paste0(original_res_back$Node1, '______', original_res_back$Sign, '______', original_res_back$Node2)
  original_res_back <- original_res_back[,c(5,4)]
  
  shuffled_res_back$edgeIDs <- paste0(shuffled_res_back$Node1, '______', shuffled_res_back$Sign, '______', shuffled_res_back$Node2)
  shuffled_res_back <- shuffled_res_back[,c(5,4)]
  
  common_PKN <- merge(common_PKN, original_res_back, by = 'edgeIDs', all.x = T)
  common_PKN$Weight <- ifelse(is.na(common_PKN$Weight), 0, common_PKN$Weight)
  
  common_PKN <- merge(common_PKN, shuffled_res_back, by = 'edgeIDs', all.x = T)
  common_PKN$Weight.y <- ifelse(is.na(common_PKN$Weight.y), 0, common_PKN$Weight.y)
  
  common_PKN$Weight.x <- as.numeric(common_PKN$Weight.x)
  common_PKN$Weight.y <- as.numeric(common_PKN$Weight.y)
  
  common_PKN$delta <- abs(common_PKN$Weight.x - common_PKN$Weight.y)
  common_PKN$percent <- gsub('/.*','',x)
  
  common_PKN$delta <- ifelse(common_PKN$edgeIDs %in% shuffled_PKN$edgeIDs, common_PKN$delta, 100)
  
  return(common_PKN)
}, original_PKN = original_PKN, original_res_back = original_res_back)

full_backward_delta_df <- as.data.frame(do.call(rbind, delta_backward_df_list))

sub_full_backward_delta_df <- full_backward_delta_df[full_backward_delta_df$percent %in% c(2,10,20,30,40,50),]

ggplot(sub_full_backward_delta_df, aes(x = delta, fill = percent)) + 
  geom_density(alpha = 0.33) +
  theme_minimal()