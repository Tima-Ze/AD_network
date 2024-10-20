require(data.table)
library(dplyr, quietly = T)

args=(commandArgs(TRUE))
file <- args[1]
input_dir=args[2]
output_dir=args[3]
print(c(file, input_dir, output_dir))

if(!dir.exists(file.path(output_dir))) {
  dir.create(file.path(output_dir),recursive = TRUE)}

for (i in paste0(file, '.txt')) {
  fread(paste0(input_dir, i), check.names = F, header = T) %>% 
    filter(abs(wTO)>=0.5) %>%
    fwrite(., paste0(output_dir, gsub('1','',file),'.txt'), sep="\t", quote = F)
  rm(i)
}
