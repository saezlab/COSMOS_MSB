library(readr)
library(igraph)
library(visNetwork)
library(stringr)

setwd("~/Dropbox/COSMOS_MSB/results/COSMOS_result/")

load("COSMOS_res_session.RData")

source("../../scripts/COSMOS_functions.R")
source("../../scripts/piano_function.R")

metab_input <- as.data.frame(
  read_csv("../../data/metab_input_COSMOS.csv"))

signaling_input <- as.data.frame(read_csv("../../../data/signaling_input_COSMOS.csv"))

source("../../scripts/format_cosmos_res_function.R")

metab_to_pubchem <- as.data.frame(read_csv("../../support/metab_to_pubchem.csv"))

cosmos_forward <- format_COSMOS_res(cosmos_res = CARNIVAL_Result_rerun,
                                    metab_mapping = metab_to_pubchem,
                                    measured_nodes = c(names(metab_input),names(signaling_input_carnival)),
                                    omnipath_ptm = read_rds("../../support/omnipath_ptm.rds"), gene_mapping = "else")

cosmos_backward <- format_COSMOS_res(cosmos_res = CARNIVAL_Result_2,
                                     metab_mapping = metab_to_pubchem,
                                     measured_nodes = c(names(metab_input),names(signaling_input_carnival)),
                                     omnipath_ptm = read_rds("../../support/omnipath_ptm.rds"), gene_mapping = "else")

full_sif <- as.data.frame(unique(rbind(cosmos_forward[[1]],cosmos_backward[[1]])))
full_att <- as.data.frame(unique(rbind(cosmos_forward[[2]][,-6],cosmos_backward[[2]][,-6])))
full_att <- full_att[full_att$AvgAct != 0,]
full_att <- unique(full_att[,c(1,6,7,8)])
full_att[duplicated(full_att$Nodes),]

full_att$Nodes <- gsub(",","_",full_att$Nodes)
full_sif$Node1 <- gsub(",","_",full_sif$Node1)
full_sif$Node2 <- gsub(",","_",full_sif$Node2)
full_sif$Weight <- as.numeric(full_sif$Weight)

full_sif <- full_sif %>% 
  group_by(Node1, Node2, Sign) %>% 
  summarise(Weight= sum(Weight))

full_sif[full_sif$Weight > 100, "Weight"] <- 100

full_att <- full_att[full_att$Nodes %in% full_sif$Node1 | full_att$Nodes %in% full_sif$Node2,]

full_att[full_att$Nodes == "CREBBP","type"] <- "TF" #misslabeled
full_att[full_att$Nodes == "KMT2A","measured"] <- "1" #misslabeled
full_att[full_att$Nodes == "GNAI1","type"] <- "Kinase"

display_node_neighboorhood(c("Adenosine_c_"),sif = full_sif, att = full_att, n = 8)

hallmarks <- gmt_to_df("../../support/h.all.v7.1.symbols.gmt_20200709.gmt")

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
