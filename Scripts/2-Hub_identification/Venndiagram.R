library(scales)
library(RColorBrewer)

#We want to make a venn diagram with the proportion of overlap of hubs between the AD and control networks

#set a directory to save the output
output_dir <- 'Results/centralities/plot/'
if(!dir.exists(output_dir)){
  dir.create(path=output_dir, recursive = T)
}

#Read core hubs
ad <- read.csv('Results/centralities/CONS_AD_coreness_hubs.csv')
ctrl <- read.csv('Results/centralities/CONS_control_coreness_hubs.csv')

#Make and save venn diagram
library(ggvenn)
data <-list(control_core=ctrl$ensembl_gene_id, AD_Core=ad$ensembl_gene_id)

ggvenn(
  data, 
  fill_color = c('#8dd3c7', '#fb8072'),
  stroke_size = 2, set_name_size =8, show_percentage = T, digits = 0,
  stroke_color = 'white',stroke_alpha = 5,text_size = 7,auto_scale = T, text_color = '#152421')

ggsave(paste0(output_dir, 'venndiagram_hubs.svg'), device = 'svg', dpi = 300, width = 10, height = 8)
