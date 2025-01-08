require(data.table)

output_dir="../Results/Raw_wTO/"

#Pass arguments including N (number of bootstrap), file (input file), and save (wTO network) to the script.
args=(commandArgs(TRUE))
N <- args[1]
file <- args[2]
#Checkpoint
print(N)
print(file)

##This calls two functions that are modified wTO source Rscripts, which are modified to accelerate the process
source("./wTO_Functions_snakemake.R")

input = read.delim(file) %>%
  as.data.frame()

wto = wTOFast_edit(Data = input, n = as.numeric(N))

wto %>% fwrite(paste0(output_dir, file), sep = "\t", quote = F)

message(paste("Saved calculations to\n", file))

