library(dplyr)
library(data.table)


args=(commandArgs(TRUE))
input <- args[1]
output <- args[2]

#Checkpoint
print(input)
print(output)


if(!dir.exists(file.path(output_dir))) {
  dir.create(file.path(output_dir),recursive = TRUE)}

# Make a empty data frame to store the details of each network
table_Raw <- data.frame()

for (i in dir(path = input, pattern = "*.txt", full.names = T)) {
  Net <- fread(i, header = T,check.names = F) 
  details_Raw <- data.frame(
    name=basename(i),
    Links_Num=nrow(Net),
    Node_Num=unique(c(Net$Node.1, Net$Node.2)) %>% length(),
    positive_links_Prec= Net %>% filter(wTO>0) %>% nrow()*100/nrow(Net),
    negative_links_Prec= Net %>% filter(wTO<0) %>% nrow()*100/nrow(Net),
    range_Positive_wTO=paste(Net %>% filter(wTO>0) %>%dplyr::select(wTO) %>%  summary(), collapse = ", "),
    range_Negative_wTO=paste(Net %>% filter(wTO<0) %>%dplyr::select(wTO)%>% summary(), collapse = ", ")
  )
  table_Raw=rbind(table_Raw, details_Raw)
  fwrite(table_Raw, paste0(output), sep="\t", quote = F)
}
