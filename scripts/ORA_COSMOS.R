library(readr)
library(biomaRt)
library(piano)
library(ggplot2)
library("org.Hs.eg.db")

setwd("~/Dropbox/COSMOS_MSB/")
source("scripts/COSMOS_functions.R")
source("scripts/piano_function.R")

#Load inputs
meta_network <- as.data.frame(
  read_csv("support/metaPKN_filtered.csv"))

prots <- unique(c(meta_network$source,meta_network$target))
prots <- prots[!grepl("XMetab",prots)]
prots <- gsub("^X","",prots)
prots <- gsub("Gene[0-9]+__","",prots)
prots <- gsub("_reverse","",prots)
prots <- gsub("EXCHANGE.*","",prots)
prots <- unique(prots)
prots <- unlist(sapply(prots, function(x){unlist(strsplit(x,split = "_"))}))

gene_mapping <- "else"

if(gene_mapping == "ensembl")
{
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  
  G_list <- getBM(filters = "entrezgene_id",
                  attributes = c('hgnc_symbol','entrezgene_id', "description"),
                  values = prots, mart = ensembl)
  gene_mapping <- G_list[,1]
  names(gene_mapping) <- G_list[,2]
  
} else
{
  entrezid <- prots
  gene_mapping <- mapIds(org.Hs.eg.db, entrezid, 'SYMBOL', 'ENTREZID')
  gene_mapping <- unlist(gene_mapping)
  gene_mapping <- gene_mapping[!is.na(gene_mapping)]
}

background <- sapply(prots, function(x, translation_dictionary){
  return(translation_dictionary[x])
},translation_dictionary = gene_mapping)
names(background) <- 1:length(background)

background <- unique(background)

hallmarks <- gmt_to_df("support/h.all.v7.1.symbols.gmt_20200709.gmt")

cosmos_att <- as.data.frame(
  read_csv("results/subnet_combined_att.csv"))

att <- cosmos_att

cosmos_att <- cosmos_att[cosmos_att$Nodes %in% background,]
cosmos_att <- unique(cosmos_att$Nodes)


hallmarks_ora <- runGSAhyper(genes = cosmos_att, universe = background, gsc = loadGSC(hallmarks))
hallmarks_ora <- as.data.frame(hallmarks_ora$resTab)
hallmarks_ora$pathway <- row.names(hallmarks_ora)

hallmarks_ora <- hallmarks_ora[order(hallmarks_ora$`Adjusted p-value`, decreasing = F),]

top_hallmark <- hallmarks_ora[1:20,c(7,1)] 
top_hallmark <- top_hallmark[order(top_hallmark$`p-value`,decreasing = T),]
top_hallmark$pathway <- gsub("HALLMARK","",top_hallmark$pathway)
top_hallmark$pathway <- gsub("_"," ",top_hallmark$pathway)
top_hallmark$`p-value` <- -log10(top_hallmark$`p-value`)
top_hallmark$pathway <- factor(top_hallmark$pathway, levels = top_hallmark$pathway)
names(top_hallmark)[2] <- "-log10(p-value)"


write_csv(top_hallmark,'results/COSMOS_result/top_hallmarks.csv')

ggplot(top_hallmark, aes(x = pathway, y = `-log10(p-value)`, fill = `-log10(p-value)`)) + geom_bar(stat = "identity") + 
  coord_flip() + 
  theme_minimal() + 
  xlab("-log10(p-value)") +
  ggtitle("COSMOS hallmark ORA") +
  scale_fill_gradient(low="grey", high="darkred")

View(att[att$Nodes %in% hallmarks[hallmarks$term == "HALLMARK_INTERFERON_GAMMA_RESPONSE","gene"],c("Nodes","Activity")])

cgp <- gmt_to_df("support/c2.cgp.v7.1.symbols.gmt")
cgp <- cgp[cgp$gene %in% background,]

cgp_ora <- runGSAhyper(genes = cosmos_att, universe = background, gsc = loadGSC(cgp))
cgp_ora <- as.data.frame(cgp_ora$resTab)
cgp_ora$pathway <- row.names(cgp_ora)

cp <- gmt_to_df("support/c2.cp.v7.1.symbols.gmt")
cp <- cp[cp$gene %in% background,]

cp_ora <- runGSAhyper(genes = cosmos_att, universe = background, gsc = loadGSC(cp))
cp_ora <- as.data.frame(cp_ora$resTab)
cp_ora$pathway <- row.names(cp_ora)
