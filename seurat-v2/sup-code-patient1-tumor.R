### R script utilizing Seurat package (v2.0.1, Butler et al., 2018) 
# to generate tSNE plots for Tumor and PBMC data for discovery patient.

library(Seurat)
library(stringr)

##-- Tumor # 21,861 genes and 7,431 cells
# Tumor_AcquiredResistance             Tumor_Before 
# 5188                                    2243 

# Load the data & Normalization according to Zheng et al., 2017
raw_dataTumor <- read.csv('Output_2018-03-12/GSE117988_raw.expMatrix_Tumor.csv.gz', 
                          header = TRUE, row.names = 1)
dim(raw_dataTumor) # 7,431 cells and 21,861 genes - already filtered
dataTumor <- log2(1 + sweep(raw_dataTumor, 2, median(colSums(raw_dataTumor))/colSums(raw_dataTumor), '*')) # Normalization
cellTypes <- sapply(colnames(dataTumor), function(x) ExtractField(x, 2, '[.]'))
cellTypes <-ifelse(cellTypes == '1', 'Tumor_Before', 'Tumor_AcquiredResistance')
table(cellTypes)
# Create Seurat object
tumor <- CreateSeuratObject(raw.data = dataTumor, min.cells = 1, min.genes = 0, project = '10x_Tumor') # already normalized
tumor # 21,861 genes and 7,431 cells

# Add meta.data (nUMI and cellTypes)
tumor <- AddMetaData(object = tumor, metadata = apply(raw_dataTumor, 2, sum), col.name = 'nUMI_raw')
tumor <- AddMetaData(object = tumor, metadata = cellTypes, col.name = 'cellTypes')

# Cluster tumor cells
tumor <- ScaleData(object = tumor, vars.to.regress = c('nUMI_raw'), model.use = 'linear', use.umi = FALSE)
tumor <- FindVariableGenes(object = tumor, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
tumor <- RunPCA(object = tumor, pc.genes = tumor@var.genes)
tumor <- RunTSNE(object = tumor, dims.use = 1:10, perplexity = 25)
TSNEPlot(tumor, group.by = 'cellTypes', colors.use = c('#EF8A62', '#67A9CF'))



