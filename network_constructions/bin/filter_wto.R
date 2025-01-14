#!/usr/bin/env Rscript

require(data.table)
library(dplyr, quietly = T)


args <- (commandArgs(TRUE))
file <- args[1]
output_dir <- args[2]
wkdir <- args[3]

# set the working directory
setwd(wkdir)
getwd()
# check arguments
print(c(file, output_dir))

# make output directory if it doesn't exist
if (!dir.exists(file.path(output_dir))) {
  dir.create(file.path(output_dir), recursive = TRUE)
}

# filter links with abs(wTO) >= 0.5
read.delim(file, check.names = F, header = T) %>%
  as.data.frame() %>%
  dplyr::filter(abs(wTO) >= 0.5) %>%
  fwrite(paste0(output_dir, basename(file)), sep = "\t", quote = F)
