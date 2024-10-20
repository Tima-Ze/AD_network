
require(data.table)
input_dir="../Data/"
output_dir="../Results/Raw_wTO/"
args=(commandArgs(TRUE))
N <- args[1]
file <- args[2]
save <- args[3]
print(N)
print(save)
print(file)

# R CMD BATCH --vanilla '--args  N=1000 save="./wTO.csv" file="input.csv"' wTO_run.out &

source("./wTO_Functions_snakemake.R")

input = read.delim(paste0(input_dir, file)) %>%
  as.data.frame()

wto = wTOFast_edit(Data = input, n = as.numeric(N))

wto %>% fwrite(paste0(output_dir, save), sep = "\t", quote = F)

message(paste("Saved calculations to\n", save))

