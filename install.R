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


remove.packages('Seurat')
pkgs = c( 'mixtools', 'lars', 'dtw', 'doSNOW', 'hdf5r' ) 
pkgs=c("ggplot2", "cowplot", "ROCR", "mixtools", "lars", "ica", "tsne", "Rtsne", "fpc", "ape", "pbapply", "igraph", "RANN", "dplyr", "RColorBrewer", "irlba", "reshape2", "gplots", "dtw", "SDMTools", "plotly", "Hmisc", "httr", "tidyr", "ggridges", "metap", "lmtest", "fitdistrplus", "png", "doSNOW", "reticulate", "foreach", "hdf5r", "RcppEigen", "RcppProgress")
#pkgs=c('jackstraw','slingshot')
BiocManager::install(pkgs,ask = F,update = F)
packageurl <- "https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_2.3.4.tar.gz"
packageurl
install.packages(packageurl, repos=NULL, type="source")
library(Seurat)

BiocManager::install('monocle')

pkgs=c('curl',openssl,hdf5r,httr,plotly)
BiocManager::install(pkgs,ask = F,update = F)

