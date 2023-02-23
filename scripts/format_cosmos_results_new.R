cleanup_protein_node <- function(node)
{
  prefixe <- str_extract(node, "Gene[0-9]+__")
  prefixe <- gsub("Gene","Enzyme",prefixe)
  
  suffixe_direction <- str_extract(node, "_reverse")
  suffixe_exchange <- str_extract(node, "EXCHANGE[0-9]+")
  
  node <- gsub("Gene[0-9]+__","",node)
  node <- gsub("_reverse","",node)
  node <- gsub("EXCHANGE[0-9]+","",node)
  
  node <- str_split(node,'_')[[1]]
  
  return(list("prefixe" = prefixe, "node" = node,"suffixe_direction" = suffixe_direction, "suffixe_exchange" = suffixe_exchange))
}

cleanup_metab_node <- function(node)
{
  prefixe <- 'Metab__'
  
  suffixe <- str_extract(node,'___[a-z]____')
  
  node <- gsub('Metab__','',node)
  node <- gsub('___[a-z]____','',node)
  
  return(list("prefixe" = prefixe, "node" = node,"suffixe_compartment" = suffixe))
}

translate_column <- function(sif_column,
                             metab_mapping,
                             gene_mapping)
{
  sif_column <- sapply(sif_column, function(x,
                                            metab_mapping,
                                            gene_mapping)
    {

    node <- x
    
    type_node <- ifelse(grepl("Metab",node),'metabolite','protein')
    
    if(type_node == 'protein')
    {
      node <- cleanup_protein_node(node)
      
      node$node <- unlist(sapply(node$node, function(x,
                                              gene_mapping)
        {
        return(ifelse(x %in% names(gene_mapping),
                      gene_mapping[x],
                      x))
      },simplify = T,
      gene_mapping = gene_mapping, USE.NAMES = F))
      
      node <- node[!is.na(node)]
      
      node$node <- paste(node$node, collapse = '_')
      
      node <- paste(node, collapse = '')
  
    } else
    {
      node <- cleanup_metab_node(node)
      
      node$node <- unlist(sapply(node$node, function(x,
                                                     metab_mapping)
      {
        return(ifelse(x %in% names(metab_mapping),
                      metab_mapping[x],
                      x))
      },simplify = T,
      metab_mapping = metab_mapping, USE.NAMES = F))
      
      node <- node[!is.na(node)]
      node <- paste(node, collapse = '')
    }
    
    
      
  }, simplify = T, 
  metab_mapping = metab_mapping,
  gene_mapping = gene_mapping, USE.NAMES = F)
  
  return(sif_column)
}

translate_sif <- function(sif,
                          metab_mapping,
                          gene_mapping)
{
  sif$Node1 <- translate_column(sif$Node1,
                                metab_mapping,
                                gene_mapping)
  
  sif$Node2 <- translate_column(sif$Node2,
                                metab_mapping,
                                gene_mapping)
  
  return(sif)
}


#' format_COSMOS_res
#'
#' formats the network with readable names
#'
#' @param cosmos_res  results of CARNIVAL run
#' @param metab_mapping a named vector with pubchem cid as names and desired metabolite names as values.
#' @param gene_mapping by default, use the 'org.Hs.eg.db' to map gene names. Can also be a named vector with entrez gene id as names and desired gene names as values.
#' @param measured_nodes vector of nodes that are measured or inputs
#' @param omnipath_ptm ptms database from OmnipathR
format_COSMOS_res <- function(cosmos_res,
                              metab_mapping,
                              gene_mapping = 'org.Hs.eg.db',
                              measured_nodes,
                              omnipath_ptm)
{
  require(dorothea)
  require('org.Hs.eg.db')
  sif <- as.data.frame(cosmos_res$weightedSIF)
  sif$Node1 <- gsub("^X","",sif$Node1)
  sif$Node2 <- gsub("^X","",sif$Node2)
  att <- as.data.frame(cosmos_res$nodesAttributes)#[,c(1,2)]
  names(att)[1] <- "Nodes"
  att$measured <- ifelse(att$Nodes %in% measured_nodes, 1, 0)
  att$Nodes <- gsub("^X","",att$Nodes)
  att$type <- ifelse(grepl("Metab",att$Nodes), "metabolite","protein")
  att <- att[abs(as.numeric(att$AvgAct)) > 0,]
  
  ########################
  prots <- unique(att$Nodes)
  prots <- prots[!(grepl("Metab",prots))]
  prots <- gsub("Gene[0-9]+__","",prots)
  prots <- gsub("_reverse","",prots)
  prots <- gsub("EXCHANGE.*","",prots)
  prots <- unique(prots)
  prots <- unlist(sapply(prots, function(x){unlist(strsplit(x,split = "_"))}))
  
  if(is.null(names(gene_mapping)))
  {
    entrezid <- prots
    gene_mapping <- mapIds(org.Hs.eg.db, entrezid, 'SYMBOL', 'ENTREZID')
    gene_mapping <- unlist(gene_mapping)
    gene_mapping <- gene_mapping[!is.na(gene_mapping)]
  } else
  {
    if(sum(names(gene_mapping) %in% prots) == 0)
    {
      print('Error : none of gene mapping identifiers are not found in the network solution. Please make sure you inputed a properlly named vector. Else use the default value for this argument : "org.Hs.eg.db"')
      return('Bad identifer mapping')
    }
  }
  
  sif <- translate_sif(sif,
                       metab_mapping,
                       gene_mapping)
  
  att$Nodes <- translate_column(att$Nodes,
                                metab_mapping,
                                gene_mapping)
  
  ########################
  omnipath_ptm <- omnipath_ptm[omnipath_ptm$modification %in% c("dephosphorylation","phosphorylation"),]
  KSN <- omnipath_ptm[,c(4,3)]
  KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_type, sep ="_")
  KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_offset, sep = "")
  KSN$sign <- ifelse(omnipath_ptm$modification == "phosphorylation", 1, -1)
  att$type <- ifelse(att$type == 'protein', 'metab_enzyme', att$type)
  
  att$type <- ifelse(att$Nodes %in% KSN$enzyme_genesymbol, "Kinase",att$type)
  dorothea <- as.data.frame(dorothea::dorothea_hs[dorothea::dorothea_hs$confidence %in% c("A","B","C","D","E"),c(3,1,4)])
  names(dorothea) <- c("target_genesymbol","source_genesymbol","sign")
  att$type <- ifelse(att$Nodes %in% dorothea$source_genesymbol, "TF",att$type)

  att$Activity <- sign(as.numeric(as.character(att$AvgAct)))

  sif <- sif[sif$Node1 != sif$Node2,]

  sif <- sif[sif$Node1 %in% att$Nodes & sif$Node2 %in% att$Nodes,]
  
  measured <- att[att$measured == 1,"Nodes"]
  att$measured <- ifelse(att$Nodes %in% measured, 1, 0)
  
  return(list(sif,att))
}