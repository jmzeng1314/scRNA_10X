### R script utilizing Seurat package (v2.0.1, Butler et al., 2018) to generate tSNE plots for Tumor and PBMC data for discovery patient.

library(Seurat)
library(stringr)

##-- Tumor
# Load the data & Normalization according to Zheng et al., 2017
raw_dataTumor <- read.csv('raw.expMatrix_Tumor.csv', header = TRUE, row.names = 1)
dim(raw_dataTumor) # 7,431 cells and 21,861 genes - already filtered
dataTumor <- log2(1 + sweep(raw_dataTumor, 2, median(colSums(raw_dataTumor))/colSums(raw_dataTumor), '*')) # Normalization
cellTypes <- sapply(colnames(dataTumor), function(x) ExtractField(x, 2, '[.]'))
cellTypes <-ifelse(cellTypes == '1', 'Tumor_Before', 'Tumor_AcquiredResistance')

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



##-- PBMC
# Load the data & Normalization according to Zheng et al., 2017
raw_dataPBMC <- read.csv('raw.expMatrix_PBMC.csv', header = TRUE, row.names = 1)
dim(raw_dataPBMC) # 12,874 cells and 17,712 genes - already filtered
dataPBMC <- log2(1 + sweep(raw_dataPBMC, 2, median(colSums(raw_dataPBMC))/colSums(raw_dataPBMC), '*')) # Normalization
timePoints <- sapply(colnames(dataPBMC), function(x) ExtractField(x, 2, '[.]'))
timePoints <-ifelse(timePoints == '1', 'PBMC_Pre', 
                    ifelse(timePoints == '2', 'PBMC_EarlyD27',
                           ifelse(timePoints == '3', 'PBMC_RespD376', 'PBMC_ARD614')))

# Create Seurat object
PBMC <- CreateSeuratObject(raw.data = dataPBMC, min.cells = 1, min.genes = 0, project = '10x_PBMC') # already normalized
PBMC # 17,712 genes and 12,874 cells

# Add meta.data (nUMI and timePoints)
PBMC <- AddMetaData(object = PBMC, metadata = apply(raw_dataPBMC, 2, sum), col.name = 'nUMI_raw')
PBMC <- AddMetaData(object = PBMC, metadata = cellTypes, col.name = 'TimePoints')

# Cluster PBMC
PBMC <- ScaleData(object = PBMC, vars.to.regress = c('nUMI_raw'), model.use = 'linear', use.umi = FALSE)
PBMC <- FindVariableGenes(object = PBMC, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
PBMC <- RunPCA(object = PBMC, pc.genes = PBMC@var.genes)
PBMC <- FindClusters(object = PBMC, reduction.type = "pca", dims.use = 1:10, resolution = 1, k.param = 35, save.SNN = TRUE) # 13 clusters
PBMC <- RunTSNE(object = PBMC, dims.use = 1:10)
TSNEPlot(PBMC, colors.use = c('green4', 'pink', '#FF7F00', 'orchid', '#99c9fb', 'dodgerblue2', 'grey30', 'yellow', 'grey60', 'grey', 'red', '#FB9A99', 'black'))



 