rm(list = ls()) # clear the environment
#load all the necessary libraries
options(warn=-1) # turn off warning message globally
suppressMessages(library(Seurat))



## 读入文章关于第一个病人的PBMC表达矩阵


start_time <- Sys.time()
raw_dataPBMC <- read.csv('Output_2018-03-12/GSE117988_raw.expMatrix_PBMC.csv.gz', header = TRUE, row.names = 1)
end_time <- Sys.time()
end_time - start_time
raw_dataPBMC[1:4,1:4]
dim(raw_dataPBMC) 

start_time <- Sys.time()
dataPBMC <- log2(1 + sweep(raw_dataPBMC, 2, median(colSums(raw_dataPBMC))/colSums(raw_dataPBMC), '*')) # Normalization
end_time <- Sys.time()
end_time - start_time


timePoints <- sapply(colnames(dataPBMC), function(x) ExtractField(x, 2, '[.]'))

timePoints <-ifelse(timePoints == '1', 'PBMC_Pre', 
                    ifelse(timePoints == '2', 'PBMC_EarlyD27',
                           ifelse(timePoints == '3', 'PBMC_RespD376', 'PBMC_ARD614')))
table(timePoints)


## 表达矩阵的质量控制
 
fivenum(apply(dataPBMC,1,function(x) sum(x>0) ))
boxplot(apply(dataPBMC,1,function(x) sum(x>0) ))
fivenum(apply(dataPBMC,2,function(x) sum(x>0) ))
hist(apply(dataPBMC,2,function(x) sum(x>0) ))



## 然后创建Seurat的对象


start_time <- Sys.time()
# Create Seurat object
PBMC <- CreateSeuratObject(raw.data = dataPBMC, 
                           min.cells = 1, min.genes = 0, project = '10x_PBMC') # already normalized
PBMC # 17,712 genes and 12,874 cells
end_time <- Sys.time()
end_time - start_time

# Add meta.data (nUMI and timePoints)
PBMC <- AddMetaData(object = PBMC, 
                    metadata = apply(raw_dataPBMC, 2, sum),
                    col.name = 'nUMI_raw')
PBMC <- AddMetaData(object = PBMC, metadata = timePoints, col.name = 'TimePoints')



## 一些质控
 
sce=PBMC
VlnPlot(object = sce, 
        features.plot = c("nGene", "nUMI"), 
        group.by = 'TimePoints', nCol = 2)
GenePlot(object = sce, gene1 = "nUMI", gene2 = "nGene")



可以看看高表达量基因是哪些


tail(sort(Matrix::rowSums(sce@raw.data)))
## 散点图可视化任意两个基因的一些属性（通常是细胞的度量）
# 这里选取两个基因。
tmp=names(sort(Matrix::rowSums(sce@raw.data),decreasing = T))
GenePlot(object = sce, gene1 = tmp[1], gene2 = tmp[2])

# 散点图可视化任意两个细胞的一些属性（通常是基因的度量）
# 这里选取两个细胞
CellPlot(sce,sce@cell.names[3],sce@cell.names[4],do.ident = FALSE)



## 最后标准聚类可视化


start_time <- Sys.time()
# Cluster PBMC
PBMC <- ScaleData(object = PBMC, vars.to.regress = c('nUMI_raw'), model.use = 'linear', use.umi = FALSE)
PBMC <- FindVariableGenes(object = PBMC, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
PBMC <- RunPCA(object = PBMC, pc.genes = PBMC@var.genes)
## 避免太多log日志被打印出来。
PBMC <- FindClusters(object = PBMC, 
                     reduction.type = "pca", 
                     dims.use = 1:10, 
                     resolution = 1, 
                     print.output = 0,
                     k.param = 35, save.SNN = TRUE) # 13 clusters
PBMC <- RunTSNE(object = PBMC, dims.use = 1:10)
TSNEPlot(PBMC, colors.use = c('green4', 'pink', '#FF7F00', 'orchid', '#99c9fb', 'dodgerblue2', 'grey30', 'yellow', 'grey60', 'grey', 'red', '#FB9A99', 'black'))
end_time <- Sys.time()
end_time - start_time
# save(PBMC,file = 'patient1.PBMC.output.Rdata')
# 这个步骤输出文件 1.75G, 遂放弃！

 
## 显示运行环境


sessionInfo()




