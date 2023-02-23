library(readr)
library(igraph)
library(visNetwork)
library(stringr)

setwd("~/Dropbox/COSMOS_MSB/results/CPTAC/")

source("../../scripts/COSMOS_functions.R")
source("../../scripts/piano_function.R")

full_sif <- as.data.frame(read_csv("subnet_combined_sif_full_newDoro_long.csv"))
full_att <- as.data.frame(read_csv("subnet_combined_att_full_newDoro_long.csv"))

hallmarks <- gmt_to_df("../../support/h.all.v7.1.symbols.gmt_20200709.gmt")
#,'MTOR','TP53','HIF1A','STAT3'
MTOR <- full_att[full_att$Nodes %in% c(hallmarks[hallmarks$term == "HALLMARK_PI3K_AKT_MTOR_SIGNALING","gene"],'MTOR','TP53','HIF1A','STAT3'),]

MTOR_sif <- full_sif[full_sif$Node1 %in% MTOR$Nodes | full_sif$Node2 %in% MTOR$Nodes,]

display_node_neighboorhood('MTOR',sif = MTOR_sif, att = MTOR, n = 10)


fatty_att <- full_att[full_att$Nodes %in% hallmarks[hallmarks$term == "HALLMARK_FATTY_ACID_METABOLISM","gene"],]

INFg <- full_att[full_att$Nodes %in% hallmarks[hallmarks$term == "HALLMARK_INTERFERON_GAMMA_RESPONSE","gene"],]

G2M <- full_att[full_att$Nodes %in% hallmarks[hallmarks$term == "HALLMARK_G2M_CHECKPOINT","gene"],]

OPHP <- full_att[full_att$Nodes %in% hallmarks[hallmarks$term == "HALLMARK_OXIDATIVE_PHOSPHORYLATION","gene"],]

#### PAPER FIGURE

node_set <- c(G2M$Nodes, INFg$Nodes,"YY1", "BCAT1","GLUL","GGT1", "ESR1","JUN","NFKB1","AKT3","SRPK2","PRKCD","Adenine_c_","Adenosine_c_",
              "L-Glutamate_c_","Reduced Glutathione_c_","L-Glutamine_c_",
              "CCNH","ADA","SMAD4","MAPK1","Inosine_c_","Hypoxanthine_c_","KMT2A","S-Adenosyl-L-Homocysteine_c_","AHCY","ADORA2B","GNAI1"
              ,"YES1","AHCYL1")

node_set <- node_set[!node_set %in% c("CDK1","MAPK14","CDC7","E2F4","MEIS1","MEIS2","CMPK2","IRF9","CIITA","IL6","TFDP1","IRF1","IRF2","HMGA1","STAT2","STAT1","UPP1","ESR1","SMAD4","AHCY","ADORA2B","GNAI1")]

test <- full_sif[full_sif$Node1 %in% node_set & full_sif$Node2 %in% node_set,]
node_set[!(node_set %in% test$Node1 | node_set %in% test$Node2)]

display_node_neighboorhood(node_set,sif = test, att = full_att, n = 1) %>% 
  visNodes(font = list(size = 25)) 

display_node_neighboorhood(c(OPHP$Nodes,'CRAT'),sif = full_sif, att = full_att, n = 2) %>% 
  visNodes(font = list(size = 25)) 

display_node_neighboorhood(c('Succinate_c_','Succinate_m_','Isovaleryl Carnitine_c_','O-Propanoylcarnitine_x_','O-Propanoylcarnitine_c_','O-Propanoylcarnitine_m_'),sif = full_sif, att = full_att, n = 3) %>% 
  visNodes(font = list(size = 25)) 
