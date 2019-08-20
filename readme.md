# 单细胞转录组10X数据处理视频教程配套代码

上游分析见：https://mp.weixin.qq.com/s/L_jpPh7TrGc-cU7EPnYhXw  

最先在生信技能树有过一些教程：

- [10X genomics单细胞数据集探索](https://mp.weixin.qq.com/s?__biz=MzI1Njk4ODE0MQ==&mid=2247483671&idx=1&sn=981cf1b1b853c03d00ba212a337793d1&scene=21#wechat_redirect)
- [10x的单细胞转录组数据就应该这样处理](https://mp.weixin.qq.com/s?__biz=MzI1Njk4ODE0MQ==&mid=2247483678&idx=1&sn=763d3ede5fd474a88fcee09f3412d7f7&scene=21#wechat_redirect)

我们技能树的学习者为这个项目专门在**单细胞天地公众号**有非常详细的教程：

- [单细胞实战(一)数据下载](https://mp.weixin.qq.com/s?__biz=MzI1Njk4ODE0MQ==&mid=2247484146&idx=1&sn=16e09b82d048eed1ff6100b22970abd5&scene=21#wechat_redirect)
- [单细胞实战(二) cell ranger使用前注意事项](https://mp.weixin.qq.com/s?__biz=MzI1Njk4ODE0MQ==&mid=2247484179&idx=1&sn=fe84f5243a6021fe6afea128e3ac273a&scene=21#wechat_redirect)
- [单细胞实战(三) Cell Ranger使用初探](https://mp.weixin.qq.com/s?__biz=MzI1Njk4ODE0MQ==&mid=2247484206&idx=1&sn=edeebbdd092f79361aee87e9ce086d80&scene=21#wechat_redirect)
- [单细胞实战(四) Cell Ranger流程概览](https://mp.weixin.qq.com/s?__biz=MzI1Njk4ODE0MQ==&mid=2247484355&idx=1&sn=7860fe0c46073a55d2d3700822c3103b&scene=21#wechat_redirect)
- [单细胞实战(五) 理解cellranger count的结果](https://mp.weixin.qq.com/s?__biz=MzI1Njk4ODE0MQ==&mid=2247484402&idx=1&sn=95c2be0dc6499e4b1eb9a91d79e584d1&scene=21#wechat_redirect)

### 本视频课程内容上游分析目录是：

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

本视频课程内容下游分析内容目录是：

- **使用seurat包进行细胞分群（计算量消耗大）**
- 根据marker对细胞亚群命名
- 使用monocle包进行差异分析
- 查看指定基因表达量
- 是否一定要使用R包，以及3大R包的异同

