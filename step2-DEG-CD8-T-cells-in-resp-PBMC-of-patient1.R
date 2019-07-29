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

# 使用monocle做差异分析

rm(list = ls()) # clear the environment
#load all the necessary libraries
options(warn=-1) # turn off warning message globally
suppressMessages(library(Seurat))
load(file = 'PBMC_RespD376_for_DEG.Rdata')
count_matrix[1:4,1:4]
dim(count_matrix)
table(cluster)

## 首先构建 monocle 对象
library(monocle) 
expr_matrix <- as.matrix(count_matrix)
sample_sheet <- data.frame(cells=names(count_matrix),  
                           cellType=cluster)
rownames(sample_sheet)<- names(count_matrix)
gene_annotation <- as.data.frame(rownames(count_matrix))
rownames(gene_annotation)<- rownames(count_matrix)
colnames(gene_annotation)<- "genes"
pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)

# Create a CellDataSet from the relative expression levels
HSMM <- newCellDataSet(
  as(expr_matrix, "sparseMatrix"),
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit=0.5,
  expressionFamily=negbinomial.size()
)

HSMM <- detectGenes(HSMM, min_expr = 1)
# HSMM <- HSMM[fData(HSMM)$num_cells_expressed > 5, ]
HSMM <- HSMM[fData(HSMM)$num_cells_expressed > 1, ]
HSMM
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM

cds=HSMM
# 单细胞转录组最重要的就是把细胞分群啦，这里可供选择的算法非常多，我们首先演示PCA结果。
# 并不是所有的基因都有作用，所以先进行挑选，合适的基因用来进行聚类。
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(cds) 
# plot_pc_variance_explained(cds, return_all = F) # norm_method='log'
# 其中 num_dim 参数选择基于上面的PCA图
cds <- reduceDimension(cds, max_components = 2, num_dim = 6,
                       reduction_method = 'tSNE', verbose = T)
cds <- clusterCells(cds, num_clusters = 5) 
plot_cell_clusters(cds, 1, 2, color = "cellType")
table(pData(cds)$Cluster,cluster)
plot_cell_clusters(cds, 1, 2 )


## 接着找差异基因

if(F){
  Sys.time()
  diff_test_res <- differentialGeneTest(cds,
                                        fullModelFormulaStr = "~cellType")
  Sys.time()
  # 可以看到运行耗时
  save(diff_test_res,file='PBMC_RespD376_diff_test_res.Rdata')
}
load(file = 'PBMC_RespD376_diff_test_res.Rdata')
# Select genes that are significant at an FDR < 10%
sig_genes <- subset(diff_test_res, qval < 0.1)
dim(sig_genes)
sig_genes$gene_short_name = rownames(sig_genes)
head(sig_genes[,c("gene_short_name", "pval", "qval")] )

##  最后推断发育轨迹，并不是文章步骤，只不过是细胞数量很少，计算量小，就跑一遍玩。
if(F){
  ## 首先挑选合适的基因
  # 这里选取统计学显著的差异基因列表
  ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
  cds <- setOrderingFilter(cds, ordering_genes)
  plot_ordering_genes(cds)
  cds <- reduceDimension(cds, max_components = 2,
                         method = 'DDRTree')
  
  # 然后降维
  cds <- reduceDimension(cds, max_components = 2,
                         method = 'DDRTree')
  # 降维是为了更好的展示数据。
  # 降维有很多种方法, 不同方法的最后展示的图都不太一样, 其中“DDRTree”是Monocle2使用的默认方法
  
  # 接着对细胞进行排序
  cds <- orderCells(cds)
  
  ## 最后两个可视化函数 
  plot_cell_trajectory(cds, color_by = "Cluster")  
  plot_cell_trajectory(cds, color_by = "cellType")  
  # 可以很明显看到细胞的发育轨迹
  
  ## 可视化差异基因的表达矩阵。
  
  ## 这里可以展现marker基因在发育轨迹推断的效果，本例子随便 选取了6个差异表达基因。
  plot_genes_in_pseudotime(cds[head(sig_genes$gene_short_name),], 
                           color_by = "Cluster")
  
}


load(file = 'PBMC_RespD376_for_DEG.Rdata')
count_matrix[1:4,1:4]
dim(count_matrix)
table(cluster)
htmapGenes=c(
  'GAPDH','CD52','TRAC','IL32','ACTB','ACTG1','COTL1',
  'GZMA','GZMB','GZMH','GNLY'
)
htmapGenes %in% rownames(sig_genes)
library(pheatmap)
dat=count_matrix[htmapGenes,]
pheatmap(dat)
n=t(scale(t(dat)))
n[n>2]=2 #限定上限，使表达量大于2的等于2
n[n< -2]= -2 #限定下限，使表达量小于-2的等于-2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
ac=data.frame(group=cluster)
rownames(ac)=colnames(n)
pheatmap(n,annotation_col = ac,
         show_colnames =F,show_rownames = T)
n[n< -1]= -1 #限定下限，使表达量小于-2的等于-2
n[1:4,1:4] 
pheatmap(n,annotation_col = ac,
         show_colnames =F,show_rownames = T)

