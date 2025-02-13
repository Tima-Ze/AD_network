#!/usr/bin/env Rscript

require(data.table)
library(dplyr)

# Pass arguments including N (number of bootstrap), file (input file), and save (wTO network) to the script.
args <- (commandArgs(TRUE))
N <- args[1] # number of bootstrap
file <- args[2]
output_dir <- args[3]
wkdir <- args[4]

# set the working directory
setwd(wkdir)
getwd()

# make output directory if it doesn't exist
if (!dir.exists(file.path(output_dir))) {
  dir.create(file.path(output_dir), recursive = TRUE)
}

## This calls two functions that are modified wTO source Rscripts, which are modified to accelerate the process
source("bin/wTO_Functions.R")

input <- read.delim(file) %>%
  as.data.frame()

wto <- wTOFast_edit(Data = input, n = as.numeric(N))

wto %>% fwrite(paste0(output_dir, basename(file)), sep = "\t", quote = F)

message(paste("Saved calculations to\n", basename(file)))
