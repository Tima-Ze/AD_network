library(dplyr)
library(ggplot2)

#We want to make a bar plot with the proportion of the different types of genes in the hubs

#set a directory to save the output
output_dir <- 'Results/centralities/plot/'
if(!dir.exists(output_dir)){
  dir.create(path=output_dir, recursive = T)
}

#Read hubs and calculate the proportion of each type of gene
ad <- read.csv('Results/centralities/CONS_AD_coreness_hubs.csv')%>% 
  group_by(gene_biotype) %>%
  summarise(count = n()) %>%
  mutate(prec = round(count / sum(count,  na.rm = TRUE) * 100)) %>% 
  mutate(Network='AD')


ctrl <- read.csv('Results/centralities/CONS_control_coreness_hubs.csv')%>% 
  group_by(gene_biotype) %>%
  summarise(count = n()) %>%
  mutate(prec = round(count / sum(count,  na.rm = TRUE) * 100))%>% 
  mutate(Network='Control')

#make and save the bar plot
plot <- bind_rows(ad, ctrl)
ggplot(plot, aes(x =Network, y = prec, fill = gene_biotype)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(title = "Gene Type Distribution within Hubs", x = "Hub Group", y = "Percentage", fill = "Gene Type") +
  theme_minimal()+
  scale_fill_manual(values = c('#e28743','#154c79'))

ggsave(plot = last_plot(), filename = paste0(output_dir,'hubs_gene_type.svg'),device = 'svg',
       dpi = 300, height = 4, width = 6)

