#remotes::install_github("alberto-valdeolivas/RandomWalkRestartMH") 
library(RandomWalkRestartMH)
library(igraph)
library(biomaRt)
library(doParallel)
library(foreach)
library(dplyr)
library(tidyr)


#set a directory to save the output
output_dir <- 'Results/random_walk/leave_out/'

if(!dir.exists(file.path(output_dir))) {
  dir.create(file.path(output_dir),recursive = TRUE)}

#list clusters' names
cluster=sub('.txt', '', list.files(path = 'Results/temp/net_clusters/', full.names = F, recursive = F))

#Run RWR for each GO categories of each cluster while leave out one gene at a time
core <- makeCluster(detectCores()[1]-3)
registerDoParallel(core)
foreach(i=cluster, 
        .packages = c('dplyr', 'igraph','biomaRt', 'RandomWalkRestartMH', 'ggplot2')) %dopar%{
          df=read.delim(list.files(path = 'Results/temp/net_clusters/',
                                   pattern = paste0(i, ".txt"),
                                   recursive = F, full.names = T)) %>%
            dplyr::select(!wTO)
          #make igraph object from clusters
          g=graph_from_data_frame(df, directed = F)
          
          #Read cluster nodes and replace NA (Nodes without enriche GOs) with 'No enriched GO'
          #and separate GO categories, one GO category per row
          node <-read.delim(list.files(path = 'Results/temp/cluster_genes/',
                                       pattern = paste0(i, '_nodes.txt'),
                                       recursive = F, full.names = T)) %>% 
            mutate_at('GO_category', ~replace(., is.na(.), 'No enriched GO'))%>%
            tidyr::separate_longer_delim(GO_category, ',')
          #Run RWR per GO category
          for(go in unique(node$GO_category)[!grepl("No enriched GO", unique(node$GO_category))]){
            SeedGene=node %>%
              filter(GO_category==go) %>% 
              pull(ensembl_gene_id)
            
            leaveOUT_df <- data.frame()
            for (n in 1:length(SeedGene)) {
              SeedNodes=SeedGene[-n]
              net_MultiplexObject <- create.multiplex(list(PPI=g))
              AdjMatrix <- compute.adjacency.matrix(net_MultiplexObject)
              AdjMatrixNorm <- normalize.multiplex.adjacency(AdjMatrix)
              RWR <- Random.Walk.Restart.Multiplex(AdjMatrixNorm,
                                                   net_MultiplexObject,SeedNodes)
              RWR <- RWR$RWRM_Results %>% 
                dplyr::rename(ensembl_gene_id=NodeNames)%>% 
                dplyr::filter(ensembl_gene_id==SeedGene[n])
              
              leaveOUT_df <- rbind(leaveOUT_df, RWR )
            }
            
            if(!dir.exists(file.path(output_dir))){
              dir.create(file.path(output_dir), recursive = T)
            }
            
            write.table(leaveOUT_df, paste0(output_dir, i, '_', go, '_leaveOUT.txt'), quote = F, sep = '\t')
            
          }
          
        }
stopCluster(core)

