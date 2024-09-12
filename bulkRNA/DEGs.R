# !/usr/bin/env Rscript

# @File       :DEGs
# @Time       :2022/8/27 00:35
# @Author     :ZhouBowen
# @Product    :DataSpell
# @Version    :R 4.1.1
# @Project    :idog
# @Description:
# @Usage      :

rm(list = ls())

#install.packages("pheatmap")
#install.packages("RColorBrewer")
#install.packages("ggplot2")
#install.packages("cowplot")
#install.packages("hexbin")
#install.packages("ggrepel")
#install.packages("ggsignif")
#BiocManager::install("DESeq2")
#BiocManager::install("vsn")

######################################
# 输入参数
project <- '/Users/zhoubw/Downloads/DEGs'
outFolder <- 'wolf_dog'
condition <- 'wolf_dog_PFC.condition.tsv'
expMatrix <- 'wolf_dog.list.count.tsv'
pValue <- 0.05
fcValue <- 2
minCounts <- 10
#####################################

# 加载包
library(DESeq2)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(vsn)
library(ggplot2)
library(cowplot)
library(hexbin)
library(stringr)
library(reshape2)
library(ggsignif)
library(dplyr)
library(tidyr)
library(ggsci)

# 设置工作目录
setwd(paste0(project))
getwd()

if (file.exists(outFolder)){
	print(paste0('output directory ', outFolder,' already exists!'))
}else{
	dir.create(outFolder)
	print(paste0('The output directory is ',outFolder,'.'))
}

in_df <- read.table(expMatrix, sep="\t", row.names=1, header=T)
#去除都是0的行
in_df <- in_df[which(rowSums(in_df) > 0),]
#去除ensembl ID的版本号
rownames(in_df)<-factor(unlist(lapply(as.character(rownames(in_df)),function(x){strsplit(x, "\\.")[[1]][1]})))
##更改列名，删除_之后的内容（即”_htseq“）
# colnames(in_df) <- gsub("\\_htseq", "", colnames(in_df))

#数据长宽转换
condition_df <- read.table(condition,header = T,fill=T,na.strings = "",sep="\t")
# condition_df <- melt(condition_df, measure.vars = c('case','control'),variable.name = 'condition',value.name ='sample')

# 判断并修正df的sample列名
correct_IDs_func <-function(df,colname){
  # 检查colname是否在列名中
  if(!colname %in% colnames(df)) {
    stop(paste0("Column name", colname, "does not exist"))
  }
  # 修正名字使其符合规范
  df[[colname]] <- gsub("-", ".", df[[colname]])
  df[[colname]] <- ifelse(grepl("^\\d", df[[colname]]), paste0("X", df[[colname]]), df[[colname]])
  # 去除空值和冗余数据
  df[df == ""] <- NA
  df <- na.omit(df) %>% distinct()
  print(head(df,3))
  return(df)
}
condition_df <- correct_IDs_func(condition_df,'sample')

#转化为因子
condition_df$condition <- factor(condition_df$condition)

#提取目标分析样本生成新矩阵
cts <- in_df[, condition_df$sample, drop = FALSE]
head(cts,3)

# 原始数据count的分布可视化
cts_long <- log2(cts+1) %>%
	pivot_longer(cols = everything(), names_to = "sample", values_to = "count")

rawCount_boxplot <- ggplot(cts_long, aes(x = sample, y = count)) +
	geom_boxplot(aes(fill = sample), color = "black") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	labs(x = "sample", y = "log2(count+1)") +
	ggtitle("rawCount") +
	theme(legend.position = "none")
ggsave(paste0(outFolder,"/rawCounts",".pdf"),rawCount_boxplot,width=8,height=4)

#构建DESeq2对象dds
dds <- DESeqDataSetFromMatrix(countData = round(cts),
							 colData = condition_df,
							 design = ~ condition)
dds

#数据过滤（默认>=10）
dim(dds)
keep <- rowSums(counts(dds)) >= minCounts
dds <- dds[keep,]
dim(dds)

#指定因子水平（对照组在前）(可选)
#dds$condition <- factor(dds$condition, levels = c("MRQ","GRD", "DRD"))
#指定对照组（可选）
dds$condition <- relevel(dds$condition, ref ="control")
dds$condition

#使用DESeq()函数进行差异分析流程；
dds <- DESeq(dds)
#使用results()函数提取分析结果；
#treated vs untreated表示log2(treated/untreated)，untreated为对照；
res <- results(dds)

#使用生成结果指定比较组的方法：
res <- results(dds, contrast=c("condition","case","control"))
res
#也可改变顺序（引起log2FC列正负号变化而已）；
#res <- results(dds, contrast=c("condition","untreated","treated"))
#如果是多个分组，如何指定特定几个比较组？
# results(dds, contrast=c("condition","treated1","untreated"))
# results(dds, contrast=c("condition","treated2","untreated"))

#name:the name of the individual effect (coefficient) 用于连续型变量；
#res <- results(dds, name="condition_treated_vs_untreated")

#并行运算；
#library("BiocParallel")
#register(MulticoreParam(4))
#dds <- DESeq(dds，parallel = TRUE)

#关于p-values和adjusted p-values；
#根据p value，对结果进行升序排列:
resOrdered <- res[order(res$pvalue),]
resOrdered

#统计adjusted p-values < 0.05 （默认）的基因数；
sum(res$padj < pValue, na.rm=TRUE)
#默认情况下，参数alpha （adjusted p value cutoff）为0.1，当然也可以自定义例如，设为0.05；
#res005 <- results(dds, alpha=0.05)
#summary(res005)
#sum(res05$padj < 0.05, na.rm=TRUE)


#导出差异分析结果数据表格
resOrdered <- resOrdered[complete.cases(resOrdered), ]
resOrdered <- resOrdered[resOrdered$pvalue != 0, ]
write.csv(as.data.frame(resOrdered),
          file=paste0(outFolder,"/diff_genes.csv"), quote = FALSE)


#如果只导出差异基因，可使用subset函数；
resSig <- subset(resOrdered, padj < pValue)
write.csv(as.data.frame(resSig),
          file=paste0(outFolder,"/diff_genes_padj",pValue,".csv"),quote = FALSE)

# 分析结果可视化导出
# 计算异常值
pdf(paste0(outFolder,"/outliers.pdf"), width=6,height=6)
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
dev.off()

# 独立过滤结果
print(paste0("alpha: ",metadata(res)$alpha))
print(metadata(res)$filterThreshold)

if (nrow(condition_df) <= 50) {
  #数据标准化，VST(variance stabilizing transformations)用于大样本；
  #vsd <- vst(dds, blind=FALSE)
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  rld <- rlog(dds, blind=FALSE)
  
  #使用常规log2(n + 1)的标准化方法；
  ntd <- normTransform(dds)
  
  head(assay(vsd), 3)
  head(assay(rld), 3)
  head(assay(ntd), 3)
  
  #导出标准化后的数据（可用于GSEA等下游分析）；
  write.csv(as.data.frame(assay(vsd)),
  		  file=paste0(outFolder,"/count_transformation_vst.csv"))
  write.csv(as.data.frame(assay(rld)),
  		  file=paste0(outFolder,"/count_transformation_rlog.csv"))
  write.csv(as.data.frame(assay(ntd)),
  		  file=paste0(outFolder,"/count_transformation_normTransform.csv"))
  
  #标准化数据count的分布可视化
  vsd_norm_long <- as.data.frame(assay(vsd)) %>%
  	pivot_longer(cols = everything(), names_to = "sample", values_to = "count")
  rld_norm_long <- as.data.frame(assay(rld)) %>%
  	pivot_longer(cols = everything(), names_to = "sample", values_to = "count")
  ntd_norm_long <- as.data.frame(assay(ntd)) %>%
  	pivot_longer(cols = everything(), names_to = "sample", values_to = "count")
  
  vsdCount_boxplot <- ggplot(vsd_norm_long, aes(x = sample, y = count)) +
  	geom_boxplot(aes(fill = sample), color = "black") +
  	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  	labs(x = "sample", y = "normalized count") +
  	ggtitle("vsdCount") +
  	theme(legend.position = "none")
  
  rldCount_boxplot <- ggplot(rld_norm_long, aes(x = sample, y = count)) +
  	geom_boxplot(aes(fill = sample), color = "black") +
  	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  	labs(x = "sample", y = "normalized count") +
  	ggtitle("rldCount") +
  	theme(legend.position = "none")
  
  ntdCount_boxplot <- ggplot(ntd_norm_long, aes(x = sample, y = count)) +
  	geom_boxplot(aes(fill = sample), color = "black") +
  	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  	labs(x = "sample", y = "normalized count") +
  	ggtitle("ntdCount") +
  	theme(legend.position = "none")
  
  pdf(paste0(outFolder,"/normalizationCounts.pdf"), width=13,height=4)
  plot_grid(vsdCount_boxplot,
  					rldCount_boxplot,
  					ntdCount_boxplot,
  					ncol = 3)
  dev.off()

} else {
  
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  
  write.csv(as.data.frame(assay(vsd)),
            file=paste0(outFolder,"/count_transformation_vst.csv"))
  
  vsd_norm_long <- as.data.frame(assay(vsd)) %>%
    pivot_longer(cols = everything(), names_to = "sample", values_to = "count")
  
  vsdCount_boxplot <- ggplot(vsd_norm_long, aes(x = sample, y = count)) +
    geom_boxplot(aes(fill = sample), color = "black") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "sample", y = "normalized count") +
    ggtitle("vsdCount") +
    theme(legend.position = "none")
  
  pdf(paste0(outFolder,"/normalizationCounts.pdf"), width=13,height=4)
  plot_grid(vsdCount_boxplot)
  dev.off()
  
} 
  
#绘制样本距离聚类热图；
#计算样本距离矩阵；
sampleDists <- dist(t(assay(vsd)))

#由向量转成矩阵；
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$condition
colnames(sampleDistMatrix) <- vsd$sample

write.csv(as.data.frame(sampleDistMatrix),
          file=paste0(outFolder,"/sample_clusters.csv"))

