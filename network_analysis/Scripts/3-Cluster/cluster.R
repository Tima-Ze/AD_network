library(igraph)
library(lsa)
library(dplyr)
library(biomaRt)
set.seed(10403)

#Make a biomart function to get gene names and gene biotype from Ensembl database
ensembl <- ensembl <- useEnsembl(biomart = "genes",
                                 dataset = "hsapiens_gene_ensembl", mirror = "asia")

biomart <- function(gene){
  getBM(attributes= c("hgnc_symbol",'gene_biotype',
                      'ensembl_gene_id'),
        filters = "ensembl_gene_id", values = gene, mart = ensembl, useCache = FALSE)
}

#set a directory to save the output
outdir <- 'Results/clusters/'
if(!dir.exists(file.path(outdir))){
  dir.create(file.path(outdir), recursive = T)
}

#Read gene list from AD pathway of KEGG database that previously downloaded
kegg=read.delim('Data/kegg/AD_kegg_genes.txt')

#Read networks and find clusters using igraph package and Louvain algorithm

for (f in list.files('Data/Networks', pattern = '.txt', recursive = F, full.names = T)) {
  df=data.table::fread(f) %>% dplyr::rename(weight=CN)
  
  g=graph_from_data_frame(df, directed = F) #Make and igraph object from the networks
  cluster=cluster_louvain(g, weights = NA, resolution = 1)
  
  for (i in as.numeric(names(sizes(cluster)[sizes(cluster) >= 20]))) {
    V(g)$cluster=membership(cluster)
    tmp= get.data.frame(g, what= "vertices")
    genes=unique(subset(tmp, cluster==i)[,'name'])
    
    g_sub=induced_subgraph(g, genes) %>% 
      igraph::simplify(edge.attr.comb = 'random')
    
    sub_graph <- igraph::as_data_frame(g_sub, what = 'edges')%>% 
      dplyr::rename(Node.1=from, Node.2=to, wTO=weight)
    
    
    nodes=biomart(genes)
    nodes=nodes %>% mutate(kegg=ifelse(ensembl_gene_id%in%kegg$ensembl_gene_id, "kegg", "No_kegg"))
    
    if(any(grepl('AD', f))){
      name="AD_cluster"
    }else{
      name="ctrl_cluster"
    }
    
    write.csv(nodes, file=paste0(outdir,name, i, '_nodes.csv'),row.names = F)
    write.table(sub_graph, paste0(outdir,name, i, '.txt'), sep='\t', quote = F, row.names = F)
  }
}

