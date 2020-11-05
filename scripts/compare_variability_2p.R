library(readr)
library(ggplot2)

setwd("~/Dropbox/COSMOS_MSB/results/variability/")

edge_df_frequency <- as.data.frame(read_csv("edge_df_frequency.csv"))

sif_2p <- readRDS('../shuffler/2/full_sif_2p.Rds')

sif_2p$edgeIDs <- paste0(sif_2p$Node1,'______',sif_2p$Sign,'______',sif_2p$Node2)
sif_2p$is_present <- TRUE

comparison_edges <- merge(edge_df_frequency,sif_2p[,c(5,6)], by = 'edgeIDs', all.x = T) 
comparison_edges$is_present <- ifelse(!is.na(comparison_edges$is_present), 'present', 'absent')
comparison_edges <- comparison_edges[order(comparison_edges$frequency_percent, decreasing = F),]
comparison_edges$index <- factor(1:length(comparison_edges[,1]), levels = 1:length(comparison_edges[,1]))

ggplot(comparison_edges, aes(x = index, y = frequency_percent, fill = is_present)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  theme_minimal() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_x_discrete(breaks = NULL)

ggplot(comparison_edges, aes(x = frequency_percent, color = is_present)) +
  geom_density(alpha = 0.5, size = 4) +
  theme_minimal() + 
  # scale_x_continuous(breaks = NULL) +
  # scale_y_continuous(breaks = NULL)
