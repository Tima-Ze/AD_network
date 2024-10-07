library(biomaRt)
library(doParallel)
library(foreach)

#set a directory to save the output
output_dir <- 'Results/centralities/'

if(!dir.exists(file.path(output_dir))) {
  dir.create(file.path(output_dir),recursive = TRUE)}

#First, we need to add gene names and gene biotype from Ensembl database to the nodes of the networks.
#For that we make a function to be used in the foreach loop.

ensembl <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
biomart <- function(gene){
  getBM(attributes= c("hgnc_symbol",'gene_biotype',
                      'ensembl_gene_id'),
        filters = "ensembl_gene_id", values = gene, mart = ensembl, useCache = FALSE)
}

##calculate centralities (Degree, normalized degree, betweenness, normalized betweenness, Eigenvector,
#Coreness) in parallel

core <- makeCluster(detectCores()[1]-3)
registerDoParallel(core)

net=c('CONS_AD', 'CONS_control')

foreach(i=net, .packages = c('dplyr', 'igraph', 'biomaRt')) %dopar%{
  df <- data.table::fread(list.files('Data/Networks/',
                                     pattern = paste0(i, ".txt"), recursive = F, full.names = T))%>%
    rename(weights=CN) %>% 
    dplyr::select(1:3)
  
  g=graph_from_data_frame(df, directed = F) %>% 
    set_edge_attr('weights', value = df$weights) #make graph objects from networks
  Degree=degree(g, loops = F, normalized =F)
  Degree_norm <- degree(g, loops = F, normalized = T)
  Eigenvector <- eigen_centrality(g, directed = F, scale = T)$vector
  Closeness <- closeness(g,normalized = T)
  Betweenness<- betweenness(g, directed = F, normalized = F)
  Betweenness_norm<- betweenness(g, directed = F, normalized = T)
  Coreness <- coreness(g, mode = "all")
  centralities <- cbind(Degree, Degree_norm, Eigenvector,  Closeness, Betweenness, Betweenness_norm, Coreness) %>% 
    data.frame() %>% 
    tibble::rownames_to_column(var='ensembl_gene_id')
  
  gene=biomart(centralities$ensembl_gene_id)
  centralities <- merge(centralities, gene, by = 'ensembl_gene_id') %>%
    dplyr::select(ensembl_gene_id, hgnc_symbol,gene_biotype, Coreness,
                  Degree, Degree_norm, Betweenness, Betweenness_norm, Eigenvector,          
                  Closeness)
  write.csv(centralities, file = paste0(output_dir, i,'_centralities.csv'), row.names = F)
}
stopCluster(core)




##Pick top hubs based on maximum coreness scores
for(i in list.files(output_dir, full.names = T, recursive = F, pattern = '.csv')){
  read.csv(i) %>% 
    slice_max(Coreness, with_ties = TRUE) %>% 
    write.csv(.,paste0(output_dir,sub('_centralities','_coreness_hubs', basename(i))), row.names = F)
}

