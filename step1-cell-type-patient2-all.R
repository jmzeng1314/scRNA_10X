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
load('~/Documents/10x/patient2.seurat.output.Rdata')
end_time <- Sys.time()
end_time - start_time 
# 16G内存的MAC pro, 30.50766 secs
# 32G内存的MAC 台式机
TSNEPlot(seurat, group.by = 'cellTypes', colors.use = c('#EF8A62', '#67A9CF'))
ggsave(filename = 'TSNEPlot_patient1_all.pdf')
TSNEPlot(seurat,group.by = "ident",pt.shape ='cellTypes') 


allGenes = row.names(seurat@raw.data)
markerGenes <- c(
  "NCAM1",
  "KRT20",
  "MKI67",
  "CD3D",
  "CD3E",
  "IL7R",
  "CD8A",
  "CCR7",
  "SELL",
  "FOXP3",
  "CTLA4",
  "NKG7",
  "GZMA",
  "MS4A1",
  "CD79A",
  "CD14",
  "FCGR3A",
  "CST3",
  "FCER1A",
  "PPBP",
  "HBA1")
markerGenes %in% allGenes


# Visualize canonical marker genes as violin plots.
pdf('patient2_marker_VlnPlot.pdf', width=10, height=15)
VlnPlot(object = seurat, features.plot = markerGenes , 
        use.raw = TRUE, y.log = TRUE) 
dev.off()

# Visualize canonical marker genes on the sctransform embedding.
pdf('patient2_marker_FeaturePlot.pdf', width=10, height=15)
FeaturePlot(object = seurat, 
            features.plot =markerGenes, 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
dev.off()


## 根据这两幅图对细胞进行命名。
table(seurat@ident)
a=read.table('celltype-patient2-all.txt')
labers=a[match(as.numeric(seurat@ident),a[,1]),2]
labers=as.character(labers)
table(labers)
dim(seurat@raw.data)
seurat <- AddMetaData(object = seurat, 
                    metadata = labers, 
                    col.name = 'labers')
tmp=seurat@meta.data
table(seurat@meta.data$labers)
# 不知道为什么AddMetaData 失效。 
seurat@meta.data$labers=labers
tmp=seurat@meta.data
TSNEPlot(seurat, group.by = 'labers',
         colors.use =  colP,
         do.label = T)
# 需要修改颜色的顺序
table(labers)
head(labers)
labers=as.factor(labers)
colP=colP[match(levels(labers),a[,2])]
head(labers)
seurat@meta.data$labers=labers 

TSNEPlot(seurat,group.by = "ident",pt.shape ='cellTypes') 
TSNEPlot(seurat, group.by = 'labers',
         #cols.use = c("grey", "blue"), 
         do.label = T)
table(labers)
table(seurat@ident)
# 可以看到，这里的命名是错的，留给大家做思考题。
ggsave(filename = 'TSNEPlot_patient1_all_wrong.pdf')

 
## 接下来按照细胞类型进行拆分绘图
# PBMC,cluster: 1,2,4,6,9,13,10
# Tumor,cluster:0,3,5,7,8,11,12,14
cellTypes = seurat@meta.data$cellTypes
table(cellTypes)
seurat_PBMC = SubsetData(seurat,cellTypes =='PBMC')
TSNEPlot(seurat_PBMC, 
        do.label = T)
ggsave('tSNE_patient2-PBMC.pdf')

seurat_Tumor= SubsetData(seurat,cellTypes =='Tumor')
TSNEPlot(seurat_Tumor, 
         do.label = T)   
ggsave('tSNE_patient2-Tumor.pdf')




