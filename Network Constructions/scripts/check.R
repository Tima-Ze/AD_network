library(dplyr)
library(data.table)

input_dir='../cns'
output_dir='../'
print(c(file,input_dir, output_dir))

if(!dir.exists(file.path(output_dir))) {
  dir.create(file.path(output_dir),recursive = TRUE)}

table_Raw <- data.frame()

for (i in dir(path = input_dir, pattern = "*.txt", full.names = T)) {
  Net <- fread(i, header = T,check.names = F) %>% 
    dplyr::rename(wTO=CN)
  details_Raw <- data.frame(
    name=basename(i),
    Network_size=nrow(Net),
    Network_nodeNum=unique(c(Net$Node.1, Net$Node.2)) %>% length(),
    positive_links_Prec= Net %>% filter(wTO>0) %>% nrow()*100/nrow(Net),
    negative_links_Prec= Net %>% filter(wTO<0) %>% nrow()*100/nrow(Net),
    wTO_range_Positive_wTO=paste(Net %>% filter(wTO>0) %>%dplyr::select(wTO) %>%  summary(), collapse = ", "),
    wTO_range_Negative_wTO=paste(Net %>% filter(wTO<0) %>%dplyr::select(wTO)%>% summary(), collapse = ", ")
  )
  table_Raw=rbind(table_Raw, details_Raw)
  fwrite(table_Raw, paste0(output_dir,'cns_Networks_check.txt'), sep="\t", quote = F)
}
