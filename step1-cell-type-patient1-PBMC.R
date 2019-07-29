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

#  
rm(list = ls()) # clear the environment
#load all the necessary libraries
options(warn=-1) # turn off warning message globally
suppressMessages(library(Seurat))
# 首先加载前面使用Seurat包进行细胞分群的结果
start_time <- Sys.time()
load('~/Documents/10x/patient1.PBMC.output.Rdata')
end_time <- Sys.time()
end_time - start_time
# 16G内存的MAC pro
# Time difference of 28.80561 secs
colP<-c('green4', 
        'pink', 
        '#FF7F00', 
        'orchid', 
        '#99c9fb', 
        'dodgerblue2', 
        'grey30', 
        'yellow', 
        'grey60', 
        'grey', 
        'red', 
        '#FB9A99', 
        'black'
)
TSNEPlot(PBMC, 
         colors.use =  colP,
         do.label = T)
ggsave(filename = 'TSNEPlot_patient1_PBMC.pdf')
table(PBMC@meta.data$TimePoints)
# 可以看到作者对4个数据集的合并做的非常棒！
TSNEPlot(PBMC,group.by = "TimePoints")
table(PBMC@meta.data$TimePoints,PBMC@ident)

## 然后根据文章，可视化那些marker基因
allGenes = row.names(PBMC@raw.data)
markerGenes <- c(
  "CD3D",
  "CD3E",
  "TRAC",
  "IL7R",
  "GZMA",
  "FCGR3A",
  "CD14",
  "MS4A1",
  "FCER1A" 
)
markerGenes %in% allGenes
# Visualize canonical marker genes as violin plots.
pdf('patient1_pBMC_marker_VlnPlot.pdf', width=10, height=15)
VlnPlot(object = PBMC, features.plot = markerGenes , 
        use.raw = TRUE, y.log = TRUE) 
dev.off()

# Visualize canonical marker genes on the sctransform embedding.
pdf('patient1_pBMC_marker_FeaturePlot.pdf', width=10, height=15)
FeaturePlot(object = PBMC, 
            features.plot =markerGenes, 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
dev.off()

## 根据这两幅图对细胞进行命名。
head(PBMC@ident)
head(as.character(PBMC@ident))
head(as.numeric(as.character(PBMC@ident)))
head(as.numeric(PBMC@ident))

tmp=PBMC@meta.data
a=read.table('celltype-patient1-PBMC.txt')
labers=a[match(as.numeric(as.character(PBMC@ident)),a[,1]),2]
table(labers)
table(PBMC@ident)

dim(PBMC@raw.data)
PBMC <- AddMetaData(object = PBMC, 
                    metadata = labers, 
                    col.name = 'labers')
# 不知道为什么AddMetaData 失效。
tmp=PBMC@meta.data
table(PBMC@meta.data$labers)
PBMC@meta.data$labers=labers
tmp=PBMC@meta.data
TSNEPlot(PBMC, group.by = 'labers',
         colors.use =  colP,
         do.label = T)
# 需要修改颜色的顺序
table(labers)
head(labers)
labers=as.factor(labers)
colP=colP[match(levels(labers),a[,2])]
head(labers)
PBMC@meta.data$labers=labers 
TSNEPlot(PBMC, group.by = 'labers',
         colors.use =  colP,
         do.label = T)
ggsave(filename = 'TSNEPlot_patient1_PBMC_new.pdf')


## 接下来按照时间点进行拆分绘图

TimePoints = PBMC@meta.data$TimePoints
table(TimePoints)
PBMC_ARD614 = SubsetData(PBMC,TimePoints =='PBMC_ARD614')
TSNEPlot(PBMC_ARD614, 
         colors.use = c('green4', 'pink', '#FF7F00', 'orchid', '#99c9fb', 'dodgerblue2', 'grey30', 'yellow', 'grey60', 'grey', 'red', '#FB9A99', 'black'),
         do.label = T)
ggsave('PBMC_ARD614_PBMC_tSNE.pdf')

PBMC_EarlyD27    = SubsetData(PBMC,TimePoints =='PBMC_EarlyD27')
TSNEPlot(PBMC_EarlyD27, 
         colors.use = c('green4', 'pink', '#FF7F00', 'orchid', '#99c9fb', 'dodgerblue2', 'grey30', 'yellow', 'grey60', 'grey', 'red', '#FB9A99', 'black'),
         do.label = T)   
ggsave('PBMC_EarlyD27_PBMC_tSNE.pdf')

PBMC_Pre  = SubsetData(PBMC,TimePoints =='PBMC_Pre')
TSNEPlot(PBMC_Pre, 
         colors.use = c('green4', 'pink', '#FF7F00', 'orchid', '#99c9fb', 'dodgerblue2', 'grey30', 'yellow', 'grey60', 'grey', 'red', '#FB9A99', 'black'),
         do.label = T)
ggsave('PBMC_Pre_PBMC_tSNE.pdf')

PBMC_RespD376 = SubsetData(PBMC,TimePoints =='PBMC_RespD376')
TSNEPlot(PBMC_RespD376, 
         colors.use = c('green4', 'pink', '#FF7F00', 'orchid', '#99c9fb', 'dodgerblue2', 'grey30', 'yellow', 'grey60', 'grey', 'red', '#FB9A99', 'black'),
         do.label = T)
ggsave('PBMC_RespD376_PBMC_tSNE.pdf')

table(TimePoints)
table(PBMC_Pre@ident)
table(PBMC_EarlyD27@ident)
table(PBMC_RespD376@ident)
table(PBMC_ARD614@ident)

#  save cluster 4 and 10 for PBMC_RespD376
PBMC_RespD376@ident
PBMC_RespD376_for_DEG = SubsetData(PBMC_RespD376,
                                   PBMC_RespD376@ident %in% c(4,10))
# save(PBMC_RespD376_for_DEG,file = 'PBMC_RespD376_for_DEG.Rdata')
count_matrix=PBMC_RespD376_for_DEG@data
count_matrix[1:4,1:4]
cluster=PBMC_RespD376_for_DEG@ident
table(cluster)
save(count_matrix,cluster,
     file = 'PBMC_RespD376_for_DEG.Rdata')












