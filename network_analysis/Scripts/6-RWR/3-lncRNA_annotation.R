library(dplyr)
library(patchwork)

#We want to annotate lncRNAs to the GO categories based on their Rwr scores. We will use the top 10% of the scores as a cutoff to select the lncRNAs that are most likely to be associated with the GO categories.

#Make a directory to store the results
outdir <- 'Results/random_walk/'
if(!dir.exists(path=outdir)){
  dir.create(file.path(outdir), recursive = T)
}

#List names of the clusters
cluster=sub('.txt', '', list.files(path = 'Results/temp/net_clusters/', full.names = F, recursive = F))

#For each GO category of each cluster, read the Rwr scores for each lncRNA and select the top 10% based on the scores
for (i in cluster) {
  list=list.files(paste0(outdir, i), full.names = T, recursive = F, pattern = '.txt')
  all_scores <- lapply(list, function(file) {
    read.delim(file) %>%
      mutate(GO=sub('.txt', '', basename(file))) #Add a column with the GO category
  })
  df <- do.call(rbind, all_scores) #merge the scores
  sel_lnc=data.frame()
  for(go in unique(df$GO)){
    sel_lnc <- df %>% 
      filter(GO==go) %>% 
      group_by(GO) %>%
      filter(Score>=quantile(Score, 0.90)) %>% #select the top 10% of the scores
      rbind(.,sel_lnc)
    #Write the selected lncRNAs
    write.table(sel_lnc, paste0(outdir, i, '/', 'selected_lnc.txt'), sep = '\t', quote = F, row.names = F)
  }
  #We'd like to have an over view of overlap of lncRNAs annotated to same GO category
  plot <- sel_lnc %>% 
    group_by(ensembl_gene_id) %>% 
    summarise(n=n())
  png(filename = paste0(outdir, i, '/', 'Selected_lnc_overlap.png'), bg = 'white', width = 8,
      height    = 6,
      units     = "in",
      res       = 300,
      pointsize = 8)
  plot(density(plot$n),
       main=paste0('Selected lncRNAs overlap in',' ', i),
       xlab='Number of GOs')
  dev.off()
}
