library(ggplot2)

setwd("~/Dropbox/COSMOS_MSB/results/")

originalVar_edge_jc_mat_FORWARD <- readRDS("variability/originalVar_edge_jc_mat_FORWARD.rds")
originalVar_edge_jc_mat_BACKWARD <- readRDS("variability/originalVar_edge_jc_mat_BACKWARD.rds")

originalVar_edge_jc_mat_FORWARD_80 <- readRDS("variability/originalVar_edge_jc_mat_FORWARD_80.rds")
originalVar_edge_jc_mat_BACKWARD_80 <- readRDS("variability/originalVar_edge_jc_mat_BACKWARD_80.rds")


jc_FORWARD <- unlist(c(originalVar_edge_jc_mat_FORWARD))
jc_FORWARD <- jc_FORWARD[!is.na(jc_FORWARD)]
jc_FORWARD <- as.data.frame(jc_FORWARD)
FORWARD_mean <- mean(jc_FORWARD$jc_FORWARD)

ggplot(jc_FORWARD, aes(x = jc_FORWARD)) + geom_density(aes(fill = 'blue')) + geom_vline(xintercept = FORWARD_mean)+ theme_minimal()

jc_BACKWARD <- unlist(c(originalVar_edge_jc_mat_BACKWARD))
jc_BACKWARD <- jc_BACKWARD[!is.na(jc_BACKWARD)]
jc_BACKWARD <- as.data.frame(jc_BACKWARD)
BACKWARD_mean <- mean(jc_BACKWARD$jc_BACKWARD)

ggplot(jc_BACKWARD, aes(x = jc_BACKWARD)) + geom_density(aes(fill = 'blue')) + geom_vline(xintercept = BACKWARD_mean)+ theme_minimal()

jc_FORWARD_80 <- unlist(c(originalVar_edge_jc_mat_FORWARD_80))
jc_FORWARD_80 <- jc_FORWARD_80[!is.na(jc_FORWARD_80)]
jc_FORWARD_80 <- as.data.frame(jc_FORWARD_80)
FORWARD_80_mean <- mean(jc_FORWARD_80$jc_FORWARD_80)

ggplot(jc_FORWARD_80, aes(x = jc_FORWARD_80)) + geom_density(aes(fill = 'blue')) + geom_vline(xintercept = FORWARD_80_mean)+ theme_minimal()

jc_BACKWARD_80 <- unlist(c(originalVar_edge_jc_mat_BACKWARD_80))
jc_BACKWARD_80 <- jc_BACKWARD_80[!is.na(jc_BACKWARD_80)]
jc_BACKWARD_80 <- as.data.frame(jc_BACKWARD_80)
BACKWARD_80_mean <- mean(jc_BACKWARD_80$jc_BACKWARD_80)

ggplot(jc_BACKWARD_80, aes(x = jc_BACKWARD_80)) + geom_density(aes(fill = 'blue')) + geom_vline(xintercept = BACKWARD_80_mean)+ theme_minimal()

names(jc_FORWARD_80)[1] <- 'jc_FORWARD'
FORWARD_combined <- as.data.frame(rbind(jc_FORWARD, jc_FORWARD_80))
FORWARD_combined$weight_threshold <- c(rep('0',length(jc_FORWARD[,1])),rep('80',length(jc_FORWARD_80[,1])))

ggplot(FORWARD_combined, aes(x = jc_FORWARD, fill = weight_threshold)) + geom_density(alpha = 0.5) + theme_minimal()


names(jc_BACKWARD_80)[1] <- 'jc_BACKWARD'
BACKWARD_combined <- as.data.frame(rbind(jc_BACKWARD, jc_BACKWARD_80))
BACKWARD_combined$weight_threshold <- c(rep('0',length(jc_BACKWARD[,1])),rep('80',length(jc_BACKWARD_80[,1])))

ggplot(BACKWARD_combined, aes(x = jc_BACKWARD, fill = weight_threshold)) + geom_density(alpha = 0.5) + theme_minimal()







twop_shuffles_edge_jc_mat_FORWARD <- readRDS("shuffler_2p/2p_shuffles_edge_jc_mat_FORWARD.rds")
twop_shuffles_edge_jc_mat_BACKWARD <- readRDS("shuffler_2p/2p_shuffles_edge_jc_mat_BACKWARD.rds")

forward_jc_df <- as.data.frame(c(unlist(c(originalVar_edge_jc_mat_FORWARD)),unlist(c(twop_shuffles_edge_jc_mat_FORWARD))))
forward_jc_df$group <- "twop_shuffles"
forward_jc_df[1:length(unlist(c(originalVar_edge_jc_mat_FORWARD))),"group"] <- "originalVar"
forward_jc_df <- forward_jc_df[complete.cases(forward_jc_df),]
names(forward_jc_df)[1] <- "value"

ggplot(forward_jc_df, aes(x = value, fill = group)) + geom_density(alpha = 0.5) + theme_minimal()

backward_jc_df <- as.data.frame(c(unlist(c(originalVar_edge_jc_mat_BACKWARD)),unlist(c(twop_shuffles_edge_jc_mat_BACKWARD))))
backward_jc_df$group <- "twop_shuffles"
backward_jc_df[1:length(unlist(c(originalVar_edge_jc_mat_BACKWARD))),"group"] <- "originalVar"
backward_jc_df <- backward_jc_df[complete.cases(backward_jc_df),]
names(backward_jc_df)[1] <- "value"

ggplot(backward_jc_df, aes(x = value, fill = group)) + geom_density(alpha = 0.5) + theme_minimal() # + geom_vline(xintercept = vertical.lines)

