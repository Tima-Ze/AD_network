library(dplyr)
library(data.table)
library(igraph)

input_dir='filter_wto/'
output_dir='cns/'

wTO_consensus=function(x){
  cns <- x%>% 
    igraph::simplify(edge.attr.comb = 'random')
  cns <- igraph::as_data_frame(cns, what = 'edges')%>% 
    rename(Node.1=from, Node.2=to)
  
  data_x = cns[, -c(1, 2)]
  sum_x = apply(data_x, 1, sum)
  div = (data_x/sum_x) * data_x
  wTO_cons = apply(div, 1, sum)
  
  cons_wto = cns %>% mutate(CONS=wTO_cons) %>%
    mutate(CN = case_when(
      if_all(starts_with('wTO_'), ~ . > 0) ~ CONS,
      if_all(starts_with('wTO_'), ~ . < 0)  ~ CONS,
      TRUE ~ 0
    )) %>% 
    mutate(CN=round(CN, digits = 2)) %>% 
    dplyr::select(Node.1, Node.2, CN) %>% 
    filter(abs(CN)>0)
  
  return(cons_wto)
}


if(!dir.exists(file.path(output_dir))) {
  dir.create(file.path(output_dir),recursive = TRUE)}

for (i in dir(input_dir, pattern = '.txt', full.names = T)) {
  fread(i, check.names = F, header = T) %>%
    select(Node.1,  Node.2,  wTO) %>%
    graph_from_data_frame(., directed = F) %>% 
    assign(gsub('_wTO.txt', '', basename(i)), ., envir = .GlobalEnv)
  rm(i)
}


#AD_TCX
intersection(AD_mayo_TCX, AD_msbbBM22,AD_msbbBM36) %>%
  wTO_consensus() %>% 
  fwrite(., paste0(output_dir,'CONS_AD_TCX.txt'), sep="\t", quote = F)

intersection(control_mayo_TCX, control_msbbBM22,control_msbbBM36) %>% 
  wTO_consensus() %>% 
  fwrite(., paste0(output_dir,'CONS_control_TCX.txt'), sep="\t", quote = F)

#AD_FCX
intersection(AD_msbbBM44, AD_msbbBM10, byname = T, keep.all.vertices = F) %>% 
  wTO_consensus() %>% 
  fwrite(., paste0(output_dir,'CONS_AD_FCX.txt'), sep="\t", quote = F)

intersection(control_msbbBM44, control_msbbBM10) %>% 
  wTO_consensus() %>% 
  fwrite(., paste0(output_dir,'CONS_control_FCX.txt'), sep="\t", quote = F)
