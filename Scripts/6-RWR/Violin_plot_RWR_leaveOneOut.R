library(ggplot2)
library(dplyr)

#To make a violin plot of the RWR scores of the GO terms in the leave-out analysis read the scores.
input <- 'Results/random_walk/leave_out/'

#set a directory to save the output
outdir <- 'Results/random_walk/plot/'
if(!dir.exists(file.path(outdir))){
  dir.create(file.path(outdir), recursive = T)
}

#list the clusters, read and merge them into one single dataframe
for (i in cluster) {
  list=list.files(input, full.names = T, recursive = F, pattern = i)
  all_scores <- lapply(list, function(file) {
    read.delim(file) %>% 
      mutate(group=sub(".*_.*_(.*)_.*\\.txt", "\\1", basename(file)))
  })
  plot <- do.call(rbind, all_scores)%>%
    mutate(Score=round(Score, digits = 6)) ##Round the scores to 6 decimal places
  options(scipen = 999) ##To avoid scientific notation
  
  #plot the violin plot
  ggplot(plot) +
    aes(x = "", y = Score, fill = group) +
    geom_violin(adjust = 1L, scale = "area") +
    theme_minimal()+
    labs(title = i, fill="GO terms")+
    theme(strip.text= element_blank())+
    scale_fill_brewer(palette = "Set3")+
    scale_y_continuous(name = "RWR scores", limits = c(min(plot$Score), max(plot$Score)))+
    facet_wrap(~group, scales = 'fixed', ncol = length(unique(plot$group)))
  ggsave(plot=last_plot(), filename = paste0(outdir, i, '_leaveOUT.png'), bg = 'white')
}
