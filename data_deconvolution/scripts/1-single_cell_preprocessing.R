library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)

##save single cell data and meta data in data/ folder before running the script. Load them here with a appropriate Seurat function

meta=read.csv('data/snRNA-seq_cell_meta.csv')
ini_data <- Read10X(data.dir = data)

#create Seurat object
subj <- CreateSeuratObject(counts = ini_data, min.cells = 3, min.features = 200)

#QC
subj <- PercentageFeatureSet(subj, pattern = "^MT-", col.name = "percent.mt")

#filter cells based on gene counts and mitochondria percentage
subj <- subset(subj, subset = nFeature_RNA > 200 & 
                 nFeature_RNA < 7500 & 
                 percent.mt < 5)
#Normalization
subj <- SCTransform(subj, vars.to.regress = "percent.mt", verbose = FALSE)

#add meta data
subj=AddMetaData(subj, meta)

#Find clusters
subj <- FindNeighbors(subj, dims = 1:10, verbose = FALSE) %>% 
  FindClusters( resolution = 0.4) %>% 
  RunUMAP( dims = 1:10) %>% 
  RunTSNE(dims = 1:10)
DimPlot(subj, reduction = "umap", label = T)

#GENE MARKERS
#Astrocytes
FeaturePlot(subj, features = c('AQP4', 'SLC1A2', 'GPC5'),label = T )
#Oligocytes
FeaturePlot(subj, features = c('MBP', 'MOBP', 'PLP1'), label = T)
#microgli
FeaturePlot(subj, features = c('CD74', 'LPAR6', 'DOCK8', 'CSF1R'), label = T)
#Endothelial
FeaturePlot(subj, features = c('FLT1', 'CLDN5', 'LAMA2'), label = T)
#OPC
FeaturePlot(subj, features = c('VCAN', 'DSCAM', 'PCDH15', 'MEGF11'),label = T)
#neuron
FeaturePlot(subj, features = c('CNTNAP2', 'SYT1', 'RBFOX1'),label = T)
#Excitory_neuron
FeaturePlot(subj, features = 'NRGN', label = T)
#Inhibitory_neuron
FeaturePlot(subj, features = 'GAD1', label = T)


#Assign cell type identity to clusters
subj=RenameIdents(object = subj, '0'= "Oligodendrocytes (Oli)", '1'="Oligodendrocytes (Oli)",
                  "2"="Oligodendrocytes (Oli)",'3'= "Oligodendrocytes (Oli)",'13'= "Oligodendrocytes (Oli)",
                  "5"="Astrocytes (Ast)","12"="Astrocytes (Ast)",
                  '10'="Microglia (Mic)",'7'="Microglia (Mic)",
                  "6"="Oligodendrocyte progenitors (Opc)",
                  "16"='Endothelial (End)', 
                  "4"='Excitatory neurons(Ex)', "11"='Excitatory neurons(Ex)',
                  "14"='Excitatory neurons(Ex)', "17"='Excitatory neurons(Ex)', 
                  "8"='Inhibitory neurons (in)', "9"='Inhibitory neurons (in)', "15"='Inhibitory neurons (in)',
                  "18"='Inhibitory neurons (in)')


DimPlot(subj, reduction = "umap", label = TRUE, pt.size = 0.5)
#save Seurat object
saveRDS(subj, file = "data/single_cell_clusters.rds")   
