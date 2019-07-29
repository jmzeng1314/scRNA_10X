rm(list = ls()) # clear the environment
#load all the necessary libraries
options(warn=-1) # turn off warning message globally
suppressMessages(library(Seurat))
library(ggfortify)

load(file = 'PBMC_RespD376_for_DEG.Rdata')
count_matrix[1:4,1:4]
dim(count_matrix)
table(cluster) 
table(cluster)
col=rainbow(2)[as.factor(cluster)]
table(col)
## only tSNE
if(T){
  choose_T_counts=count_matrix
  choose_T_counts=choose_T_counts[apply(choose_T_counts,1, sd)>0,]
  choose_T_counts=choose_T_counts[names(tail(sort(apply(choose_T_counts,1, sd)),1000)),]
  pca_dat <- prcomp(t(choose_T_counts), scale. = TRUE)
  p=autoplot(pca_dat,col=col) + theme_classic() + ggtitle('PCA plot')
  print(p)
  
  str(pca_dat)
  pca_dat$x[1:4,1:4]
  new_pca <- FactoMineR::PCA(
    t(choose_T_counts), 
    ncp = ncol(choose_T_counts), 
    graph=FALSE
  )
  set.seed(42)
  # TSNE即t-distributed Stochastic Neighbor Embedding.
  library(Rtsne)
  tsne_out <- Rtsne(pca_dat$x[,1:5], perplexity = 10,
                    pca=FALSE, 
                    max_iter=2000, 
                    verbose=TRUE) # Run TSNE
  tsnes=tsne_out$Y
  colnames(tsnes) <- c("tSNE1", "tSNE2")
  ggplot(tsnes, aes(x = tSNE1, y = tSNE2))+ geom_point(col=col)
  
}

##　for seurat 2
if(T){
  PBMC <- CreateSeuratObject(raw.data = count_matrix,  
                             min.cells = 1, min.genes = 0, 
                             project = '10x_PBMC') # already normalized
  # Add meta.data (nUMI and cluster)
  PBMC <- AddMetaData(object = PBMC, 
                      metadata = apply(count_matrix, 2, sum),
                      col.name = 'nUMI_raw')
  PBMC <- AddMetaData(object = PBMC, 
                      metadata = cluster, 
                      col.name = 'cluster')
  sce=PBMC
  VlnPlot(object = sce, 
          features.plot = c("nGene", "nUMI"), 
          group.by = 'cluster', nCol = 2)
  tmp=sce@meta.data
  
  # 需要搞清楚去除了什么 ScaleData
  PBMC <- ScaleData(object = PBMC, 
                    vars.to.regress = c('nGene',"nUMI"), 
                    model.use = 'linear', 
                    use.umi = FALSE)
  PBMC <- FindVariableGenes(object = PBMC, 
                            mean.function = ExpMean, 
                            dispersion.function = LogVMR, 
                            x.low.cutoff = 0.0125, 
                            x.high.cutoff = 3, y.cutoff = 0.5)
  length( PBMC@var.genes)
  PBMC <- RunPCA(object = PBMC, 
                 pc.genes = PBMC@var.genes)
  PBMC <- RunTSNE(object = PBMC, dims.use = 1:20, 
                  perplexity = 25)
  TSNEPlot(PBMC, group.by = 'cluster')
  TSNEPlot(PBMC )
}


# for monocle 
if(T){
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
  cds <- clusterCells(cds, num_clusters = 3) 
  plot_cell_clusters(cds, 1, 2, color = "cellType")
  table(pData(cds)$Cluster,cluster)
  plot_cell_clusters(cds, 1, 2 )
}

## for scater 

if(T){
  library(scater)
  ct=count_matrix
  pheno_data <-  data.frame(cluster)
  pheno_data$Cell=rownames(pheno_data)
# data("sc_example_counts")
# data("sc_example_cell_info") 
# ct=sc_example_counts
# pheno_data=sc_example_cell_info
  # 参考 https://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/vignette-intro.R
dim(ct)  
dim(pheno_data)
sce <- SingleCellExperiment(
    assays = list(counts = as.matrix(ct)), 
    colData = pheno_data
  )
  sce
  # 举例
  exprs(sce) <- log2(
    calculateCPM(sce ) + 1)
  sce <- runPCA(sce)
  # 这里并没有进行任何基因的挑选，就直接进行了PCA，与 seurat包不一样。
  reducedDimNames(sce)
  ## --------------
  plotReducedDim(sce, use_dimred = "PCA", 
                 colour_by = "cluster" )
  
  
}









