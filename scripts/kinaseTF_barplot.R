library(readr)
library(ggplot2)

setwd("~/Dropbox/COSMOS_MSB/")

kinase_activities <- as.data.frame(
  read_csv("data/phospho_kinase_activities.csv"))

TF_scores <- as.data.frame(
  read_csv("data/RNA_TF_scores.csv"))

kinase_activities <- kinase_activities[order(abs(kinase_activities$NES), decreasing = T),,drop = F]

TF_scores <- TF_scores[order(abs(TF_scores$NES), decreasing = T),,drop = F]

top_kinase <- kinase_activities[1:25,]
names(top_kinase)[1] <- "ID"
top_TF <- TF_scores[1:25,]

top_kinase$type <- "kinase/phosphatase"
top_TF$type <- "TF"

top_combined <- as.data.frame(rbind(top_TF,top_kinase))
top_combined <- top_combined[order(abs(top_combined$NES), decreasing = T),]
top_combined$ID <- factor(top_combined$ID, levels = top_combined$ID)

ggplot(top_combined, aes(x = ID, y = NES, fill = type)) + 
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

