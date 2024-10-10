library(dplyr)
library(purrr)
library(stringr)

#To examine the differential co-expressed links and doing RWR on clusters, we need to do some data preprationa
#We save the preprocessed data in temp folder. #Note: Keep in mind that we only keep those clusters which have at least one enriched GO term.
#These pre-processing steps are:
#1- Identifying Enriched GOs belong to each cluster
#2-We want to add GO_categories attribute to the gene tables belong to each cluster
#3-We want to add differential co-expressed link types (from CoDiNA result) to the network link tables


################################# Let's get started ###########################################

##Read GO enrichment results
GO=read.csv('Results/GO/GO_enrichmetn.csv')

################################# 1-Annotating GOs to clusters ################################
#set a directory to save the output
output_dir <- 'Results/GO/cluster_GOs'
if(!dir.exists(output_dir)){
  dir.create(path=output_dir, recursive = T)
}

for (x in c('AD', 'control')) {
  tmp <- GO%>%
    filter(stringr::str_detect(Cluster, x)) %>% 
    filter(!(GO_category%in%c('To be removed')))%>%
    group_by(Cluster) %>%
    summarise(across(where(is.numeric), ~ sum(., na.rm = TRUE)),
              across(where(is.character), ~ paste(sort(unique(.)), collapse = ","))) %>% 
    mutate(geneID = stringr::str_replace_all(geneID, ",", "/")) %>% 
    na.omit() %>% 
    assign(x,., envir = .GlobalEnv)
}

if(!dir.exists(file.path(output_dir))){
  dir.create(file.path(output_dir), recursive = T)
}


for(D in c('AD', 'control')){
  for (i in unique(get(D)$Cluster)) {
    get(D) %>% filter(Cluster==i) %>% 
      write.table(., paste0(output_dir,i, '.txt'), sep = '\t', quote = F, row.names = F)
    rm(i)}
}


####################### 2-Add GO_category attribution to genes of clusters ##########################
#set a directory to save the output
output_dir <- 'Results/temp/cluster_genes/'
if(!dir.exists(file.path(output_dir))){
  dir.create(file.path(output_dir), recursive = T)
}

# read nodes belong to each cluster
for(i in list.files(path='Results/clusters', pattern = '*.csv', full.names = F)){
  df=read.csv(paste0('Results/clusters/', i))
  cluster <- sub('_nodes.csv','',basename(i))
  
  
  tmp <- GO%>%
    filter(stringr::str_detect(Cluster, cluster)) %>% 
    group_by(GO_category) %>%
    summarise(across(where(is.character), ~ paste(sort(unique(.)), collapse = ","))) %>% 
    mutate(geneID = stringr::str_replace_all(geneID, ",", "/")) %>% 
    filter(!(GO_category==c('To be removed'))) %>%
    na.omit()
  
  if(nrow(tmp)>=1){
    genes <- tmp %>%
      dplyr::select(geneID, GO_category,Description, Cluster) %>% 
      tidyr::separate_rows(geneID, sep='/') %>% 
      distinct() %>% 
      group_by(geneID) %>%
      summarise_all(~ paste(unique(.), collapse = ","), .groups = "drop")%>% 
      mutate(GO_category=
               str_split(GO_category, ',') %>% map(sort) %>% map_chr(paste, collapse = ','))
    
    full_join(df, genes, by=c('hgnc_symbol'='geneID'))%>%
      tidyr::replace_na(list(Cluster=cluster)) %>% 
      write.table(., paste0(output_dir, sub('.csv','',basename(i)), '.txt'), sep = '\t', quote = F, row.names = F)
  }
}

####################### 3- Add link types attribute to the each cluster edge table #################
library(igraph)
library(future)
library(doParallel)

input_dir <- 'Results/clusters/'

#set a directory to save the output
output_dir <- 'Results/temp/net_clusters/'
if(!dir.exists(file.path(output_dir))){
  dir.create(file.path(output_dir), recursive = T)
}

#Read CoDiNA differential co-expressed links
codina <- data.table::fread('Results/codina/DiffNet.txt') %>% 
  dplyr::select(Node.1, Node.2, Phi, Phi_tilde)

#Read clusters' edge tables
cluster <- sub('_nodes.txt', '', list.files(path = 'Results/temp/cluster_genes', full.names = F, recursive = F))

#We use Union function from igraph to merge CoDiNA edge attributes to our clusters' edge tables
core <- makeCluster(detectCores()[1]-3)
registerDoParallel(core)

foreach(i=cluster, .packages = c('dplyr', 'igraph'))%dopar%{
  net <- read.delim(list.files(path=input_dir, pattern = paste0(i,'.txt'), full.names = T))
  nodes <- unique(c(net$Node.1, net$Node.2))
  g1 <- graph_from_data_frame(codina, directed = F)
  g2 <- graph_from_data_frame(net, directed = F)
  g_union <- union(g1,g2)
  g_sub <- subgraph(g_union, vids = nodes)
  df <- as_data_frame(g_sub, what = 'edges') %>%
    rename(Node.1=from, Node.2=to) %>% 
    na.omit()
  
  write.table(df,paste0(output_dir, i, '.txt'), sep = '\t',  quote = F, row.names = F)
}

stopCluster(core)




