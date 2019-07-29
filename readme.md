# 单细胞转录组10X数据处理上游流程

### 背景自学

最先在生信技能树有过一些教程：

- [10X genomics单细胞数据集探索](https://mp.weixin.qq.com/s/zQK9k6oC4nzZ9uOzo-8Vmw)
- [10x的单细胞转录组数据就应该这样处理](https://mp.weixin.qq.com/s/xDc5a36oqZp4-KtzfTRoag)

我们技能树的学习者为这个项目专门在**单细胞天地公众号**有非常详细的教程：

- [单细胞实战(一)数据下载](https://mp.weixin.qq.com/s/PoqArKLtIslN2cJtXZTzOg)
- [单细胞实战(二) cell ranger使用前注意事项](https://mp.weixin.qq.com/s/fP8f4HboMM7m2Nd7AIljlg)
- [单细胞实战(三) Cell Ranger使用初探](https://mp.weixin.qq.com/s/6Jqu-20HasfHen6vRUSSBQ)
- [单细胞实战(四) Cell Ranger流程概览](https://mp.weixin.qq.com/s/v2S8obShNRpeTRFQt2PrwQ)
- [单细胞实战(五) 理解cellranger count的结果](https://mp.weixin.qq.com/s/_VGFGmYBJmYm_4KLc9zamg)

### 本文数据

下面的数据及数据库存放在华为云，两个月内有效(到2019年8月31失效)

文章解读见：（2019年3月份）第12周（总第60周 ）- 10个单细胞转录组数据探索免疫治疗机理 <http://www.bio-info-trainee.com/4089.html> 

Patient 2586-4

- 介绍：<https://www.ncbi.nlm.nih.gov//geo/query/acc.cgi?acc=GSE117988>

- 存放在：/teach/paper/raw/P2586-4
- 大小：**31G**（打包后P2586-4.tar.gz为28G）

Patient 9245-3

- 介绍：<https://www.ncbi.nlm.nih.gov//geo/query/acc.cgi?acc=GSE118056>
- 存放在：/teach/paper/raw/P9245-3
- 大小： **129G**（压缩打包108G）

### 配置cellranger软件及数据库

- **biosoft 存放软件** 
  - 准备cellranger 2.0.2软件=》分析文章数据
    - <https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/2.0/>
    - 压缩文件732M，**解压后1.8G**
  - 准备cellranger 3.0.2软件=》分析测试数据
    - <https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/3.0/>
    - 压缩文件955M，**解压后2.2G**

- **database 存放参考数据**

  - 为cellranger 2.0准备
    - cd /teach/database/cranger2
    - wget <http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-1.2.0.tar.gz>
    - 11G大小，**解压后17G**
  - 为cellranger 3.0准备
    - cd /teach/database/cranger3
    - wget <http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz>
    - 11G大小，**解压后17G**

配置好的软件和数据库如下：

```
root@singlecell:/home/t1# ls /teach/database/cranger2/refdata-cellranger-GRCh38-1.2.0
README.BEFORE.MODIFYING  fasta  genes  pickle  reference.json  star  version
root@singlecell:/home/t1# ls /teach/biosoft/cellranger-2.0.2/cellranger
/teach/biosoft/cellranger-2.0.2/cellranger
root@singlecell:/home/t1# ls -lh /teach/biosoft/cellranger-2.0.2/cellranger
lrwxrwxrwx 1 t1 teach 34 Sep  8  2017 /teach/biosoft/cellranger-2.0.2/cellranger -> cellranger-cs/2.0.2/bin/cellranger

root@singlecell:/home/t1# ls -lh /teach/paper/raw/P2586-4
total 31G
-rwxr-xr-x 1 t1 teach  2.6G Jul 15 20:53 SRR7722937.sra
-rwxr-xr-x 1 t1 teach  417M Jul 15 21:12 SRR7722937_S1_L001_I1_001.fastq.gz
-rwxr-xr-x 1 t1 teach  889M Jul 15 20:11 SRR7722937_S1_L001_R1_001.fastq.gz
-rwxr-xr-x 1 t1 teach  2.4G Jul 15 20:45 SRR7722937_S1_L001_R2_001.fastq.gz
-rwxr-xr-x 1 t1 teach  3.7G Jul 15 21:05 SRR7722938.sra
-rwxr-xr-x 1 t1 teach  593M Jul 15 21:28 SRR7722938_S1_L001_I1_001.fastq.gz
-rwxr-xr-x 1 t1 teach  1.3G Jul 15 20:36 SRR7722938_S1_L001_R1_001.fastq.gz
-rwxr-xr-x 1 t1 teach  3.4G Jul 15 20:07 SRR7722938_S1_L001_R2_001.fastq.gz
-rwxr-xr-x 1 t1 teach  1.8G Jul 15 21:22 SRR7722939.sra
```

### 使用conda安装软件

```shell
# https://mirrors.tuna.tsinghua.edu.cn/help/anaconda/
# https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/ 
wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc 
## 安装好conda后需要设置镜像。
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda
conda config --set show_channel_urls yes

conda create -n 10x
conda activate 10x
conda install -y -c bioconda   sra-tools
```

成功进入 10x 软件环境，就可以使用 prefetch和fastq-dump命令。

### 下载sra文件

```shell
prefetch SRR7722942
# 自行补充成为批量脚本
```

### 转换成为fq文件

```shell
fastq-dump --gzip --split-files -A PBMC_Disc  SRR7722942.sra
```

需要修改名字 (这个非常重要，看PPT操作)

### 对10X的fq文件运行cellranger的counts流程

拿一个数据集做测试：

```shell
ref=/teach/database/cranger2/refdata-cellranger-GRCh38-1.2.0
cr=/teach/biosoft/cellranger-2.0.2/cellranger
id=SRR7722938
## 需要有3个符合规则的fq文件，下面的代码才能运行
$cr count --id=$id \
--transcriptome=$ref \
--fastqs=/teach/paper/raw/P2586-4 \
--sample=$id \
--nosecondary \
--localcores=15 \
--localmem=30
```

报告解读。

### 导入R使用seurat分析

看GitHub代码

