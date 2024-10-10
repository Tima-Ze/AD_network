library(dplyr)
library(ReactomePA)
library(enrichplot)
library(ggplot2)
library(clusterProfiler)
library(pheatmap)
library(DOSE)
library(enrichplot)
library(org.Hs.eg.db)

#set a directory to save the output
output_table <- 'Results/GO/'
output_plot <- 'Results/GO/plot/'

#we need EntrezID to perform the GO enrichment analysis, and to filter the protein-coding genes
for (i in list.files(path = 'Results/clusters', pattern = '*nodes.csv', recursive = TRUE, full.names = TRUE)) {
  df <- read.csv(i)
  nodes <- clusterProfiler::bitr(df$ensembl_gene_id, fromType='ENSEMBL', toType='ENTREZID',
                                 OrgDb='org.Hs.eg.db', drop = TRUE)
  df <- left_join(df, nodes, by=c('ensembl_gene_id'='ENSEMBL'))%>%
    distinct(ensembl_gene_id, .keep_all=T) %>% 
    filter(gene_biotype%in%'protein_coding')
  assign(sub('_nodes.csv','',basename(i)), df, envir = .GlobalEnv)
}

#As background gene, we will use the protein-coding genes from the inputted gene counts to construct the networks
for(i in list.files(path='Data/wTO_inputs/', pattern = '*.txt' ,full.names = T)){
  Bg <- character(0)
  Bg <- union(Bg, rownames(read.table(i, stringsAsFactors = F)))
  rm(i)
}
Bg=clusterProfiler::bitr(Bg, fromType='ENSEMBL', toType=c('GENETYPE', 'ENTREZID'),
                         OrgDb='org.Hs.eg.db', drop = TRUE) %>% 
  filter(GENETYPE%in%'protein-coding')

#we want to compare the gene sets across the clusters using the compareCluster function
#We need to create a list with the gene sets of each cluster
data <- list()
for (df_name in ls(pattern = ".*_cluster.*")) {
  df <- get(df_name)
  data[[df_name]] <- df$ENTREZID
}

#Go enrichment comparison across clusters
ego <- compareCluster(data, fun="enrichGO",
                      OrgDb= org.Hs.eg.db,
                      ont           = 'bp',
                      pAdjustMethod = "fdr",
                      pvalueCutoff  = 0.05,
                      readable      = TRUE,
                      univers=Bg$ENTREZID)

#We used these GO terms to make the dot plot
dotplot(ego,font.size=13, label_format=40, 
        showCategory=c(c('cytoplasmic translation', 'ribonucleoprotein complex biogenesis','post-transcriptional regulation of gene expression',
                         'cytoplasmic translation', 'mitochondrial gene expression',
                         'apoptotic mitochondrial changes','regulation of intrinsic apoptotic signaling pathway', 'apoptotic mitochondrial changes',
                         'lysosomal lumen acidification', 'establishment of protein localization to membrane', 
                         'establishment of protein localization to membrane','protein import into mitochondrial matrix', 'mitochondrial transport',
                         'generation of precursor metabolites and energy','purine-containing compound metabolic process', 'ATP metabolic process',
                         'aerobic electron transport chain', 'cellular respiration',
                         'synaptic vesicle lumen acidification','synaptic vesicle maturation', 'mitochondrial membrane organization')))

ggsave(paste0(output_plot,'GO_clusters.png'),
       plot = last_plot(), width = 7, height = 8)
ggsave(paste0(output_plot,'GO_clusters.svg'),
       plot = last_plot(), width = 7, height = 8, dpi = 300)

write.csv(ego@compareClusterResult, paste0(output_table, 'GO_enrichment.csv'), row.names = F)

#pathway analysis 

pathway <- compareCluster(data, fun="enrichPathway",organism = "human",
                          pvalueCutoff = 0.05, pAdjustMethod = "fdr", readable = T,universe= Bg$ENTREZID)

dotplot(pathway, font.size=13, label_format=40)

ggsave(paste0(output_plot,'pathway_clusters.png'),
       plot = last_plot(), width = 7, height = 9)

ggsave(paste0(output_plot,'pathway_clusters.svg'),
       plot = last_plot(), width = 7, height = 9, dpi = 300, device = 'svg')


