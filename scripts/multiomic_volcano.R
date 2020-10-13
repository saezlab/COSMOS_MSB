library(readr)
library(ggplot2)

setwd("~/Dropbox/COSMOS_MSB/")
source("scripts/multiomic_volcano.R")

ttop_RNA <- as.data.frame(
  read_csv("data/RNA_ttop_tumorvshealthy.csv"))

ttop_phospho <- as.data.frame(
  read_csv("data/phospho_ttop_tumorVsHealthy.csv"))

ttop_metab <- as.data.frame(
  read_csv("data/metab_ttop_tumour_vs_healthy.csv"))

ttop_RNA$groups <- "RNA"
ttop_phospho$groups <- "phospho"
ttop_metab$groups <- "metabs"

ttop <- as.data.frame(rbind(ttop_RNA,ttop_phospho))
ttop <- as.data.frame(rbind(ttop, ttop_metab))

ttop$my_alpha <- 0.075
ttop[ttop$groups == "phospho","my_alpha"] <- 0.1
ttop[ttop$groups == "metabs","my_alpha"] <- 1

volcano_groups(ttop, FCIndex = 2, pValIndex = 5, IDIndex = 1)
