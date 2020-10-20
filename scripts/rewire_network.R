library(BiRewire)
library(readr)

setwd("/home/ad234505/COSMOS_revisions/rewire")

dsg <- birewire.load.dsg("metaPKN.tsv")
dsg <- dsg[-1,]
dsg$V2 <- ifelse(dsg$V2 == 1, "+", "-")
dsg <- birewire.induced.bipartite(dsg)


dsg2 <- birewire.rewire.dsg(dsg, max.iter.pos = 1000, max.iter.neg = 1000)

saveRDS(dsg2, "rewired_net_1000_1000.Rds")