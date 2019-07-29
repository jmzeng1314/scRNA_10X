## 
### ---------------
###
### Create: Jianming Zeng
### Date: 2019-07-24 15:03:19
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2019-07-24  First version
###
### ---------------

# 单独检查基因的表达量情况

rm(list = ls()) # clear the environment
#load all the necessary libraries
options(warn=-1) # turn off warning message globally
suppressMessages(library(Seurat))
# 首先加载前面使用Seurat包进行细胞分群的结果
start_time <- Sys.time()
load('~/Documents/10x/patient1.tumor.output.Rdata')
end_time <- Sys.time()
end_time - start_time
# 21 秒
TSNEPlot(tumor, group.by = 'cellTypes', colors.use = c('#EF8A62', '#67A9CF'))

count_matrix=tumor@data
count_matrix[1:4,1:4]
cluster=tumor@meta.data$cellTypes
table(cluster)


allGenes = row.names(tumor@raw.data)
allGenes[grep('HLA',allGenes)]
FeaturePlot(object = tumor, 
            features.plot ='HLA-A', 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
table(count_matrix['HLA-A',]>0, cluster)
# 可以看到 HLA-A 在免疫治疗前后并没有显著差异
FeaturePlot(object = tumor, 
            features.plot ='HLA-B', 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
table(count_matrix['HLA-B',]>0, cluster)
# 可以看到 HLA-B在免疫治疗前后差异很大。


