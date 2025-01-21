#!/usr/bin/env Rscript

library(dplyr)
library(data.table)
library(igraph)

# Pass arguments 
args <- (commandArgs(TRUE))
input_dir <- args[1]
output_dir <- args[2]
wkdir <- args[2]

# set the working directory
setwd(wkdir)
getwd()

# make output directory if it doesn't exist
if (!dir.exists(file.path(output_dir))) {
  dir.create(file.path(output_dir), recursive = TRUE)
}

for (i in dir(input_dir, pattern = '.txt', full.names = T)) {
  fread(i, check.names = F, header = T) %>%
    select(Node.1,  Node.2,  wTO) %>%
    graph_from_data_frame(., directed = F) %>% 
    assign(gsub('.txt', '', basename(i)), ., envir = .GlobalEnv)
  rm(i)
}

source("bin/cns_weight_function.R")

#treat_Consensus
intersection(treat1, treat2) %>%
  wTO_consensus() %>% 
  fwrite(., paste0(output_dir,'cns_treat.txt'), sep="\t", quote = F)


#control_Consensus
intersection(control1, control2) %>% 
  wTO_consensus() %>% 
  fwrite(., paste0(output_dir,'cns_control.txt'), sep="\t", quote = F)

