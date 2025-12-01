library(CIBERSORT)
library(doParallel)
library(foreach)
load('scripts/CIBERSORT_function.R')

# Load signature matrix and bulk expression data files
signature_matrix <- read.table("data/signature_matrix.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
bulk_files <- list.files(path = "../network_constructions/data/", pattern = "*.txt", full.names = TRUE, recursive = F)

# Create output directory if it doesn't exist
dir.create("result/", showWarnings = FALSE)


# Set up parallel computing and run the CIBERSORT algorithm
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Run CIBERSORT in parallel
results_list <- foreach(file = bulk_files, .packages = c("CIBERSORT")) %dopar% {
  result <- cibersort(signature_matrix, file, perm = 1000, QN = FALSE)
  output_file <- paste0("results/", basename(file), "_cell_composition.txt")
  write.table(result, output_file, sep = "\t", quote = FALSE, col.names = NA)
  return(result)
}

# Stop parallel processing
stopCluster(cl)