require(data.table)
library(dplyr)

# Pass arguments including N (number of bootstrap), file (input file), and save (wTO network) to the script.
args <- (commandArgs(TRUE))
N <- args[1]
file <- args[2]
wkdir <- args[3]
# Check arguments
print(N)
print(file)

# set the working directory
setwd(wkdir)

# set output directory
output_dir <- "/Results/Raw_wTO/"

## This calls two functions that are modified wTO source Rscripts, which are modified to accelerate the process
source("bin/wTO_Functions.R")

input <- read.delim(paste0("Data/", file)) %>%
  as.data.frame()

wto <- wTOFast_edit(Data = input, n = as.numeric(N))

wto %>% fwrite(paste0(output_dir, file), sep = "\t", quote = F)

message(paste("Saved calculations to\n", file))
