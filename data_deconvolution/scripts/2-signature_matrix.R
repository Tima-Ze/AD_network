library(Seurat)
library(dplyr)

# Create output directory if it doesn't exist
dir.create("data/", showWarnings = FALSE)

# Load the Seurat object
seurat_obj <- readRDS("data/single_cell_clusters.rds")
control_cells <- subset(seurat_obj, subset = Diagnosis == "Control")

#The exp data is normalized using SCTransform (SCT)
# Ensure SCT is set as the active assay
DefaultAssay(seurat_obj) <- "SCT"

#Generate the signature matrix using SCT-normalized data
signature_matrix <- AverageExpression(control_cells, group.by = "ident", assays = "SCT")$SCT
signature_matrix <- signature_matrix %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Gene")

# Save the signature matrix to a text file
write.table(signature_matrix, "data/signature_matrix.txt", sep = "\t", quote = FALSE, row.names = F)
