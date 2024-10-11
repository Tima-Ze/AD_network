library(dplyr)
library(ggplot2)

#We want to add the cut-off of lncRNA annotation to the violin plot that we previously made to show the RWR distribution 
input <- 'Results/random_walk/'

#Make a output directory
output <- 'Results/random_walk/plot/cut_off/'
if(!dir.exists(path=output)){
  dir.create(file.path(output), recursive = T)
}

#List all the clusters' names
cluster=sub('.txt', '', list.files(path = 'Results/temp/net_clusters/', full.names = F, recursive = F))

#Read scores and merge them
for (i in cluster) {
  file_list <- list.files(paste0(input, i), full.names = TRUE, recursive = FALSE, pattern = '.txt')
  file_list <- file_list[!grepl('selected_lnc', file_list)]
  
  all_scores <- lapply(file_list, function(file) {
    read.delim(file) %>% 
      mutate(group = sub('.txt', '', basename(file))) 
  })
  
  plot_data <- do.call(rbind, all_scores) %>%
    mutate(Score = round(Score, digits = 6)) # Round the score to 6 decimal places
  
  # Calculate the 90th percentile for each group
  percentile_95 <- plot_data %>%
    group_by(group) %>%
    summarise(percentile_95 = quantile(Score, 0.90))
  
  #Make violin plot with a red cut-off line
  gg <- ggplot(plot_data) +
    aes(x = "", y = Score, fill = group) +
    geom_violin(adjust = 1L, scale = "area") +
    geom_hline(data = percentile_95, aes(yintercept = percentile_95), color = "red", linetype = "dashed") + # Add horizontal lines
    theme_minimal() +
    labs(title = paste('rwr distribution and cutoff for', i), fill = "GO terms") +
    theme(strip.text = element_blank()) +
    scale_fill_brewer(palette = "Set3") +
    scale_y_continuous(name = "RWR scores", limits = c(min(plot_data$Score), max(plot_data$Score))) +
    facet_wrap(~group, scales = 'fixed', ncol = length(unique(plot_data$group)))
  
  ggsave(plot = gg, filename = paste0(output, i, '_rwr_cutoff.png'), bg = 'white')
}


