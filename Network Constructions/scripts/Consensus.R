library(dplyr)
library(data.table)
source('./consensus_function.R')

#args=(commandArgs(TRUE))
input_dir='../filter_wto/'
output_dir='../cns/'
print(c(input_dir, output_dir))

if(!dir.exists(file.path(output_dir))) {
  dir.create(file.path(output_dir),recursive = TRUE)}

print( dir(input_dir, pattern = '.txt', full.names = T))
for (i in dir(input_dir, pattern = '.txt', full.names = T)) {
  fread(i, check.names = F, header = T) %>%
    select(Node.1,  Node.2,  wTO, pval) %>% 
    mutate(pval=0) %>% 
    assign(gsub('_wTO.txt', '', basename(i)), ., envir = .GlobalEnv)
  rm(i)
}
wTO.Consensus(data=list(AD_mayo_TCX, AD_msbbBM22,AD_msbbBM36)) %>% 
fwrite(., paste0(output_dir,'CONS_AD_TCX.txt'), sep="\t", quote = F)

wTO.Consensus(data=list(control_mayo_TCX, control_msbbBM22,control_msbbBM36)) %>% 
fwrite(., paste0(output_dir,'CONS_control_TCX.txt'), sep="\t", quote = F)

wTO.Consensus(data=list(AD_msbbBM44, AD_msbbBM10)) %>% 
  fwrite(., paste0(output_dir,'AD_FCX.txt'), sep="\t", quote = F)

wTO.Consensus(data=list(control_msbbBM44, control_msbbBM10)) %>% 
  fwrite(., paste0(output_dir,'control_FCX.txt'), sep="\t", quote = F)
