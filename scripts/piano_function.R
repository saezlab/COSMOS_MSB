library(piano)
library(parallel)
library(GSEABase)
library(dplyr)

runPIANO <- function(topTable, gene_to_term, nCores = 1000, IDIndex = 1, FCIndex = 2, PvalIndex = 5, TvalIndex = 4, nPerm = 10000)
{
  
  nCores <- min(c(nCores,detectCores()-1))
  
  print(paste(nCores, " cpu(s) will be used."), sep = "")
  
  topTable <- as.data.frame(topTable)
  
  gene_to_term <- as.data.frame(gene_to_term)
  gene_to_term[,1] <- toupper(gene_to_term[,1])
  gene_to_term[,2] <- toupper(gene_to_term[,2])
  #print(head(gene_to_term))
  geneSet <- loadGSC(gene_to_term)
  
  myFC <- topTable[,FCIndex]
  names(myFC) <- toupper(topTable[,IDIndex])
  
  myPval <- topTable[,PvalIndex]
  names(myPval) <- toupper(topTable[,IDIndex])
  
  myTval <- topTable[,TvalIndex]
  names(myTval) <- toupper(topTable[,IDIndex])
  
  nPerm <- nPerm-nPerm%%nCores
  
  print(paste(as.character(nPerm), " permutations will be made (so that there is an integer x such that x*nCores=nPerm)", sep = ""))
  
  ###Run the GSA
  gsaRes1 <- runGSA(myTval, gsc=geneSet, adjMethod = "fdr", geneSetStat = "mean", ncpus = nCores, nPerm = nPerm)
  #gsaRes1 <- runGSA(myTval, gsc=geneSet, adjMethod = "fdr", geneSetStat = "mean")
  gsaRes2 <- runGSA(myTval, gsc=geneSet, adjMethod = "fdr", geneSetStat = "median", ncpus = nCores, nPerm = nPerm)
  gsaRes3 <- runGSA(myTval, gsc=geneSet, adjMethod = "fdr", geneSetStat = "sum", ncpus = nCores, nPerm = nPerm)
  gsaRes4 <- runGSA(myTval, gsc=geneSet, adjMethod = "fdr", geneSetStat = "maxmean", ncpus = nCores, nPerm = nPerm)
  gsaRes5 <- runGSA(geneLevelStats = myPval, directions = myFC, gsc=geneSet, adjMethod = "fdr", geneSetStat = "fisher", ncpus = nCores, nPerm = nPerm)
  gsaRes6 <- runGSA(geneLevelStats = myPval, directions = myFC, gsc=geneSet, adjMethod = "fdr", geneSetStat = "stouffer", ncpus = nCores, nPerm = nPerm)
  gsaRes7 <- runGSA(geneLevelStats = myPval, directions = myFC, gsc=geneSet, adjMethod = "fdr", geneSetStat = "tailStrength", ncpus = nCores, nPerm = nPerm)
  gsaRes9 <- runGSA(myTval, gsc=geneSet, adjMethod = "fdr", geneSetStat = "page", ncpus = nCores, nPerm = nPerm)
  gsaRes10 <- runGSA(geneLevelStats = myPval, directions = myFC, gsc=geneSet, adjMethod = "fdr", geneSetStat = "reporter", ncpus = nCores, nPerm = nPerm)
  gsaRes11 <- runGSA(myTval, gsc=geneSet, adjMethod = "fdr", geneSetStat = "fgsea", ncpus = nCores, nPerm = nPerm)
  
  resList <- list(gsaRes1,gsaRes2,gsaRes3,gsaRes4,gsaRes5,gsaRes6,gsaRes7,gsaRes9,gsaRes10,gsaRes11)
  names(resList) <- c("mean","median","sum","maxmean","fisher", "stouffer","tailStrength","page", "reporter","fgsea")
  
  ch <- consensusHeatmap(resList,cutoff=50000,method="median", ncharLabel = 1000, cellnote = "medianPvalue", cex = 0.2, plot = FALSE) ##The results are strange
  
  consensus <- ch$pMat
  
  return(list(consensus,ch,resList))
}

gmt_to_df <- function(gmtfile, fast = T)
{
  if(fast)
  {
    genesets = GSEABase::getGmt(con = gmtfile)
    genesets = unlist(genesets)
    
    gene_to_term =plyr::ldply(genesets,function(geneset){
      temp <- geneIds(geneset)
      temp2 <- setName(geneset)
      temp3 <- as.data.frame(cbind(temp,rep(temp2,length(temp))))
      
    },.progress = plyr::progress_text())
    names(gene_to_term) <- c("gene","term")
    return(gene_to_term[complete.cases(gene_to_term),])
  }
  else
  {
    genesets = getGmt(con = gmtfile)
    genesets = unlist(genesets)
    
    gene_to_term <- data.frame(NA,NA)
    names(gene_to_term) <- c("gene","term")
    for (geneset in genesets)
    {
      temp <- geneIds(geneset)
      temp2 <- setName(geneset)
      temp3 <- as.data.frame(cbind(temp,rep(temp2,length(temp))))
      names(temp3) <- c("gene","term")
      gene_to_term <- rbind(gene_to_term,temp3)
    }
    
    return(gene_to_term[complete.cases(gene_to_term),])
  }
}

volcano_nice <- function (df, hAss = 0.05, FCIndex, pValIndex, IDIndex, vAss = NULL,
          label = FALSE, straight = FALSE, nlabels, manual_labels = NA)
{
  df <- df[complete.cases(df), ]
  names(df)[1] <- "X1"
  hAssOri <- hAss
  hAss <- -log(hAss)
  names(df) <- gsub("adj.P.Val", "FDR", names(df))
  names(df)[FCIndex] <- "logFC"
  names(df)[pValIndex] <- "adj.P.Val"
  if (max(abs(df[, FCIndex])) >= 1) {
    xlimAbs <- ceiling(max(abs(df[, FCIndex])))
    ylimAbs <- ceiling(max(abs(-log(df[, pValIndex]))))
  }
  else {
    xlimAbs <- max(abs(df[, FCIndex]))
    ylimAbs <- max(abs(-log(df[, pValIndex])))
  }
  if (is.null(vAss)) {
    vAss <- xlimAbs/10
  }
  xneg <- function(x) abs(hAss - 1 + x/(x + vAss))
  xpos <- function(x) abs(hAss - 1 + x/(x - vAss))
  test <- function(x, y, vAss) {
    if (x < -vAss) {
      if (xneg(x) < -log(y)) {
        return("1")
      }
      else {
        return("0")
      }
    }
    else {
      if (x > vAss) {
        if (xpos(x) < -log(y)) {
          return("1")
        }
        else {
          return("0")
        }
      }
      else {
        return("0")
      }
    }
  }
  if (straight) {
    df$couleur <- ifelse(abs(df$logFC) >= vAss & df$adj.P.Val <=
                           hAssOri, "1", "0")
  }
  else {
    df$couleur <- "0"
    df$couleur <- apply(df, 1, FUN = function(x) test(as.numeric(x[FCIndex]),
                                                      as.numeric(x[pValIndex]), vAss))
  }
  df <- df[order(df$adj.P.Val, decreasing = F), ]
  df$condLabel <- df[, IDIndex]
  df[df$couleur == "0", "condLabel"] <- NA
  labels_to_keep <- c(df[c(1:nlabels), "condLabel"],manual_labels)
  df[!(df$condLabel %in% labels_to_keep), "condLabel"] <- NA
  df$couleur <- ifelse(df$couleur == "1" & df$logFC < 0, "2",
                       df$couleur)
  if (label) {
    a <- ggplot(df, aes(x = logFC, y = -log(adj.P.Val), color = couleur)) +
      geom_point(alpha = 0.5) + geom_label_repel(aes(label = condLabel)) +
      stat_function(fun = xneg, xlim = c(-xlimAbs, -vAss),
                    color = "black", alpha = 0.7) + ylim(c(0, ylimAbs)) +
      xlim(c(-xlimAbs, xlimAbs)) + stat_function(fun = xpos,
                                                 xlim = c(vAss, xlimAbs), color = "black", alpha = 0.7) +
      scale_colour_manual(values = c("grey30", "red",
                                     "royalblue3")) + theme_minimal() + theme(legend.position = "none")
  }
  else {
    if (straight) {
      a <- ggplot(df, aes(x = logFC, y = -log(adj.P.Val),
                          color = couleur)) + geom_point(alpha = 0.5) +
        geom_vline(xintercept = -vAss, color = "blue") +
        geom_vline(xintercept = vAss, color = "blue") +
        ylim(c(0, ylimAbs)) + xlim(c(-xlimAbs, xlimAbs)) +
        geom_hline(yintercept = hAss, color = "red") +
        scale_colour_manual(values = c("grey30", "red",
                                       "royalblue3")) + theme_minimal() + theme(legend.position = "none")
    }
    else {
      a <- ggplot(df, aes(x = logFC, y = -log(adj.P.Val),
                          color = couleur)) + geom_point(alpha = 0.5) +
        stat_function(fun = xneg, xlim = c(-xlimAbs,
                                           -vAss), color = "black", alpha = 0.7) + ylim(c(0,
                                                                                          ylimAbs)) + xlim(c(-xlimAbs, xlimAbs)) + stat_function(fun = xpos,
                                                                                                                                                 xlim = c(vAss, xlimAbs), color = "black", alpha = 0.7) +
        scale_colour_manual(values = c("grey30", "red",
                                       "royalblue3")) + theme_minimal() + theme(legend.position = "none")
    }
  }
  return(a)
}

volcano_genesets <- function(df, gtt, consensus, outpath, mt = NULL, get_df = 0) {
  gtt <- as.data.frame(apply(gtt, 2, function(x) toupper(x)))
  names(gtt) <- c("gene","term")
  colnames(df)[1] <- "X"
  colnames(consensus)[1] <- "X"
  
  if (!is.null(mt))
  {
    mt <- apply(mt, 2, function(x) toupper(x))
    
    names(mt) <- c("X","gene")
    
    gtt <- merge(gtt,mt, all = FALSE)
    
    names(gtt) <- c("gene","term","X")
  }
  else
  {
    names(gtt) <- c("X","term")
  }
  
  df$X <- toupper(df$X)
  df <- merge(df,gtt, all = FALSE)
  
  dir.create(outpath, recursive = T, showWarnings = F)
  setwd(outpath)
  
  terms <- consensus$X
  
  for (term in terms)
  {
    subDf <- df[df$term == term,]
    
    if (length(subDf[,1]) > 4)
    {
      print(term)
      print(length(subDf[,1]))
      a <- volcano_nice(df = subDf,FCIndex = 2,pValIndex = 5,IDIndex = 1, label = T, nlabels = 30)
      ggsave(paste(term,".pdf"), a)
    }
  }
  if (get_df == 1){
    return(df)
  }
}