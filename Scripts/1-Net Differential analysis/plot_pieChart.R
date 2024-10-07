library(dplyr)
library(stringr)
library(ggplot2)

#set a directory to save the output
output_dir <- 'Results/codina/plot/'
if(!dir.exists(output_dir)){
  dir.create(path=output_dir, recursive = T)
}

#We want to make a pie chart with the proportion of the different types of links in the network

#Read CoDiNA cancluated differential network
Diffnet <- read.delim('Results/codina/DiffNet.txt')

#We want to make two pie charts, one for the AD netwrk and one for the control network
#So we need to filter the data to get the links of the AD network and the control network
#As an example common links plus AD-specific links are the links of the AD network, we do same for the control network

for (i in c('g.AD','g.control')) {
df <- Diffnet %>% filter(Phi_tilde%in%c('a',i))

#rename links
df <- df %>%
  dplyr::select(Node.1, Node.2, Phi_tilde) %>% 
  dplyr::mutate(Links_group = case_when(
    str_detect(Phi_tilde, "g") & str_detect(Phi_tilde, "AD") ~ 'AD-specific links',
    str_detect(Phi_tilde, "g") & str_detect(Phi_tilde, "Ctrl") ~ 'Control-specific links',
    str_detect(Phi_tilde, "a") ~ 'Common links',
    str_detect(Phi_tilde, "^b") ~ 'Different signs')) 

#Calculate the proportion of each type of link to be presented in the pie chart
df <-df%>%
  dplyr::select(Links_group) %>% 
  group_by(Links_group) %>% # Variable to be transformed
  count() %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))


ggplot(df, aes(x = "", y = perc, fill = Links_group)) +
  geom_col() +
  geom_text(aes(label = labels),
            size = 10,
            position = position_stack(vjust = 0.5),
            show.legend = FALSE) +
  guides(fill = guide_legend(title = "differential co-expression links")) +
  scale_fill_manual(values = c('AD-specific links'="#F4A578", 'Common links'='#8c6bb1',
                               'Control-specific links'='#8dd3c7','Different signs'='red'))  +
  coord_polar(theta = "y") + 
  theme_void()+
  theme(legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size=20))

ggsave(plot = last_plot(), paste0(output_dir, sub('g.', '',i), '_codina.svg'), device = 'svg', dpi = 300,
       width = 10, height = 8)

ggsave(plot = last_plot(), paste0(output_dir, sub('g.', '',i), '_codina.png'), device = 'png',
       width = 10, height = 8)
}
