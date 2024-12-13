require(data.table)

input_dir="../Data/"
output_dir="../Results/Raw_wTO/"

#Pass arguments including N (number of bootstrap), file (input file/gene count table), and save (resulted wTO network) to the script.
args=(commandArgs(TRUE))
N <- args[1]
file <- args[2]
save <- args[3]
#Checkpoint
print(N)
print(save)
print(file)

##This calls two functions that are modified wTO source Rscripts, which are modified to accelerate the process
source("./wTO_Functions_snakemake.R")

input = read.delim(paste0(input_dir, file)) %>%
  as.data.frame()

wto = wTOFast_edit(Data = input, n = as.numeric(N))

wto %>% fwrite(paste0(output_dir, save), sep = "\t", quote = F)

message(paste("Saved calculations to\n", save))

