library(CARNIVAL)

load("~/Dropbox/COSMOS_MSB/results/COSMOS_result/COSMOS_res_session.RData")

sif_1 <- CARNIVAL_Result_rerun$weightedSIF
sif_2 <- CARNIVAL_Result_2$weightedSIF

sif <- as.data.frame(rbind(sif_1, sif_2))
sif <- unique(sif[,-4])

metab_input_carnival <- metab_input_carnival[,names(metab_input_carnival) %in% sif$Node1 | names(metab_input_carnival) %in% sif$Node2]
signaling_input_carnival <- signaling_input_carnival[,names(signaling_input_carnival) %in% sif$Node1 | names(signaling_input_carnival) %in% sif$Node2]

metab_input_carnival_vec <- as.numeric(metab_input_carnival[1,])
names(metab_input_carnival_vec) <- names(metab_input_carnival)

signaling_input_carnival_vec <- as.numeric(signaling_input_carnival[1,])
names(signaling_input_carnival_vec) <- names(signaling_input_carnival)

for(i in 1:3)
{
  metab_input_carnival_vec <- metab_input_carnival_vec[sample(1:length(metab_input_carnival_vec), length(metab_input_carnival_vec), replace = F)]
}
metab_input_carnival_vec <- metab_input_carnival_vec[1:10]
metab_input_carnival_vec

for(i in 1:3)
{
  signaling_input_carnival_vec <- signaling_input_carnival_vec[sample(1:length(signaling_input_carnival_vec), length(signaling_input_carnival_vec), replace = F)]
}
signaling_input_carnival_vec <- signaling_input_carnival_vec[1:15]
signaling_input_carnival_vec

metab_input_carnival <- metab_input_carnival[,names(metab_input_carnival_vec)]
signaling_input_carnival <- signaling_input_carnival[,names(signaling_input_carnival)]

CARNIVAL_Result <- runCARNIVAL(inputObj = sign(signaling_input_carnival),
                               measObj = metab_input_carnival,
                               netObj = sif,
                               solverPath = "~/Documents/cplex",
                               solver = "cplex",
                               timelimit = 3600,
                               mipGAP = 0.2)

CARNIVAL_Result_2 <- runCARNIVAL(inputObj = sign(metab_input_carnival),
                               measObj = signaling_input_carnival,
                               netObj = sif,
                               solverPath = "~/Documents/cplex",
                               solver = "cplex",
                               timelimit = 3600,
                               mipGAP = 0.2)

View(CARNIVAL_Result$weightedSIF)
View(CARNIVAL_Result_2$weightedSIF)

toy_signaling_input_carnival_vec <- signaling_input_carnival_vec
toy_metab_input_carnival_vec <- metab_input_carnival_vec
toy_sif <- sif
toy_RNA <- ttop_RNA[ttop_RNA$ID %in% toy_sif$Node1 | ttop_RNA$ID %in% toy_sif$Node2,'t']
names(toy_RNA) <- ttop_RNA[ttop_RNA$ID %in% toy_sif$Node1 | ttop_RNA$ID %in% toy_sif$Node2,'ID']
names(toy_sif) <- c("source", "interaction", "target")
toy_sif$interaction <- as.numeric(toy_sif$interaction)

omnipath_ptm <- omnipath_ptm[omnipath_ptm$modification %in% c("dephosphorylation","phosphorylation"),]

save(toy_signaling_input_carnival_vec, file = '~/Dropbox/COSMOS/data/toy_signaling_input_carnival_vec.RData')
save(toy_metab_input_carnival_vec, file = '~/Dropbox/COSMOS/data/toy_metab_input_carnival_vec.RData')
save(toy_sif, file = '~/Dropbox/COSMOS/data/toy_sif.RData')
save(toy_RNA, file = '~/Dropbox/COSMOS/data/toy_RNA.RData')

