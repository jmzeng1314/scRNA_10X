# 单细胞转录组10X数据处理视频教程配套代码

上游分析见：https://mp.weixin.qq.com/s/L_jpPh7TrGc-cU7EPnYhXw  内容目录是：

- 服务器准备（16核心64G）
- 下载软件和数据库文件
- 下载文章的sra数据
- 把sra文件转换为fastq文件并介绍
- 根据cellRanger说明书对fastq改名
- 跑cellRanger的count步骤

下游分析主要是R代码，这里简要介绍。

### 5个文件夹

一一介绍如下：

-  Output_2018-03-12 文件夹是文章附件数据，作者提供的
- figures是文章里面需要重现的图表
- output是我们代码出图后集中存放
- seurat-v2和seurat-v3是使用不同版本的Seurat包对Output_2018-03-12 文件夹文章数据走标准流程得到的细胞分群结果

### 文章图表重现

主要内容是：

- **使用seurat包进行细胞分群（计算量消耗大）**
- 根据marker对细胞亚群命名
- 使用monocle包进行差异分析
- 查看指定基因表达量
- 是否一定要使用R包，以及3大R包的异同

