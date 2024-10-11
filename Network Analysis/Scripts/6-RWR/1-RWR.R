#remotes::install_github("alberto-valdeolivas/RandomWalkRestartMH") 
library(RandomWalkRestartMH)
library(dplyr)
library(igraph)
library(doParallel)
library(foreach)


#set a directory to save the output
outdir <- 'Results/random_walk/'
if(!dir.exists(path=outdir)){
  dir.create(path=outdir, recursive = T)
}


cluster=sub('.txt', '', list.files(path = 'Results/temp/net_clusters/', full.names = F, recursive = F))

core <- makeCluster(detectCores()[1]-3)
registerDoParallel(core)
foreach(i=cluster, 
        .packages = c('dplyr', 'igraph', 'RandomWalkRestartMH')) %dopar%{
          
          df=read.delim(list.files(path = 'Results/temp/net_clusters/',
                                   pattern = i,
                                   recursive = F, full.names = T)) %>%
            dplyr::select(1:2)
          
          g=graph_from_data_frame(df, directed = F) #Make a graph object from each cluster's edge table
          
          #Read the nodes file: 1-those are not annotated to GO ters are labeled as 'No enriched GO'
          #2-expand the table by separating the GO terms
          #3-For each GO term, find the protein-coding genes that are annotated to that term to be used for seed genes in RWR
          node <-read.delim(list.files(path = 'Results/temp/cluster_genes/',
                                       pattern = paste0(i, '_nodes.txt'),
                                       recursive = F, full.names = T)) %>% 
            mutate_at('GO_category', ~replace(., is.na(.), 'No enriched GO'))
          node=node %>%
            tidyr::separate_longer_delim(GO_category, ',')
          for(go in unique(node$GO_category)[!grepl("No enriched GO", unique(node$GO_category))]){
            print(go)
            SeedGene=node %>%
              filter(GO_category==go & gene_biotype=='protein_coding') %>% 
              pull(ensembl_gene_id)
            SeedGene=intersect(unique(c(df$Node.1,df$Node.2)), SeedGene)
            
            #Run RWR 
            net_MultiplexObject <- create.multiplex(list(PPI=g))
            AdjMatrix <- compute.adjacency.matrix(net_MultiplexObject)
            AdjMatrixNorm <- normalize.multiplex.adjacency(AdjMatrix)
            RWR_results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm,
                                                         net_MultiplexObject,SeedGene)
            RWR <- RWR_results$RWRM_Results %>%
              dplyr::rename(ensembl_gene_id=NodeNames)
            
            #For the RwR calculation we consider all genes in the cluster but then we filter scores for lncRNAs
            RWR <- merge(RWR, node, by = 'ensembl_gene_id', all.y = F) %>% 
              filter(gene_biotype=='lncRNA') %>% 
              arrange(Score) %>%
              select_if(~ !any(is.na(.))) %>% 
              dplyr::select(!GO_category)
            
            if(!dir.exists(file.path(paste0(outdir, i)))){
              dir.create(file.path(paste0(outdir, i)), recursive = T)
            }
            write.table(RWR, paste0(outdir, i, '/', go, '.txt'), quote = F, sep = '\t', row.names = F)
          }
        }
stopCluster(core)


