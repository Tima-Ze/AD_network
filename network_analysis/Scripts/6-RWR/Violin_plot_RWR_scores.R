library(dplyr)
library(ggplot2)

# We want to show the distribution of lncRNAs' RWR scores for each GO category in each cluster by violin plot.

#set a output directory
input <- 'Results/random_walk/'
outdir <- 'Results/random_walk/plot/RWR/'
if(!dir.exists(path=outdir)){
  dir.create(file.path(outdir), recursive = T)
}

#List clusters with enriched GO terms
cluster=sub('.txt', '', list.files(path = 'Results/temp/net_clusters/', full.names = F, recursive = F))

#Make violin plots for each cluster
for (i in cluster) {
  list=list.files(paste0(input, i), full.names = T, recursive = F, pattern = '.txt')
  list <- list[!grepl("selected_lnc\\.txt$", list)]
  all_scores <- lapply(list, function(file) {
    read.delim(file) %>% 
      mutate(group=sub('.txt', '', basename(file))) %>% 
      mutate(Score=round(Score, digits = 6))
  })
  plot <- do.call(rbind, all_scores)
  options(scipen = 999)
  ggplot(plot) +
    aes(x = "", y = Score, fill = group) +
    geom_violin(adjust = 1L, scale = "area") +
    theme_minimal()+
    labs(title = i, fill="GO terms")+
    theme(strip.text= element_blank(),
          axis.text.x = element_blank())+
    scale_fill_brewer(palette = "Set3")+
    scale_y_continuous(name = "RWR scores", limits = c(min(plot$Score), max(plot$Score)))+
    facet_wrap(~group, scales = 'fixed', ncol = length(unique(plot$group)))
  ggsave(plot=last_plot(), filename = paste0( outdir, i, 'RWR_scores.png'), bg = 'white')
}
