# !/usr/bin/env Rscript

# @File       :DEGs.R
# @Time       :2022/8/27 00:35
# @Author     :ZhouBowen
# @Product    :DataSpell
# @Version    :R 4.1.1
# @Project    :bulkRNA_pipeline
# @Description:calculate and visualise differentially expressed genes
# @Usage      :Rscript DEGs.R

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
# parameters
project <- '/PATH/DEGs'
outFolder <- 'group_id'
condition <- 'group_id.condition.tsv'
expMatrix <- 'group_id.list.count.tsv'
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

# set working directory
setwd(paste0(project))
getwd()

if (file.exists(outFolder)){
	print(paste0('output directory ', outFolder,' already exists!'))
}else{
	dir.create(outFolder)
	print(paste0('The output directory is ',outFolder,'.'))
}

in_df <- read.table(expMatrix, sep="\t", row.names=1, header=T)
# Remove rows that are all zero
in_df <- in_df[which(rowSums(in_df) > 0),]
# Remove version number of ensembl ID
rownames(in_df)<-factor(unlist(lapply(as.character(rownames(in_df)),function(x){strsplit(x, "\\.")[[1]][1]})))

# colnames(in_df) <- gsub("\\_htseq", "", colnames(in_df))

# Data Aspect Conversion
condition_df <- read.table(condition,header = T,fill=T,na.strings = "",sep="\t")
# condition_df <- melt(condition_df, measure.vars = c('case','control'),variable.name = 'condition',value.name ='sample')

# fix the sample column name for dataframe
correct_IDs_func <-function(df,colname){
  # Check if colname is in the column
  if(!colname %in% colnames(df)) {
    stop(paste0("Column name", colname, "does not exist"))
  }
  # fix the sample name
  df[[colname]] <- gsub("-", ".", df[[colname]])
  df[[colname]] <- ifelse(grepl("^\\d", df[[colname]]), paste0("X", df[[colname]]), df[[colname]])
  # remove nulls and redundant data
  df[df == ""] <- NA
  df <- na.omit(df) %>% distinct()
  print(head(df,3))
  return(df)
}
condition_df <- correct_IDs_func(condition_df,'sample')

# transformation to factors
condition_df$condition <- factor(condition_df$condition)

# extracting target analysis samples to generate new matrices
cts <- in_df[, condition_df$sample, drop = FALSE]
head(cts,3)

# visualisation of the distribution of the raw count
cts_long <- log2(cts+1) %>%
	pivot_longer(cols = everything(), names_to = "sample", values_to = "count")

rawCount_boxplot <- ggplot(cts_long, aes(x = sample, y = count)) +
	geom_boxplot(aes(fill = sample), color = "black") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	labs(x = "sample", y = "log2(count+1)") +
	ggtitle("rawCount") +
	theme(legend.position = "none")
ggsave(paste0(outFolder,"/rawCounts",".pdf"),rawCount_boxplot,width=8,height=4)

# Constructe DESeq2 object dds
dds <- DESeqDataSetFromMatrix(countData = round(cts),
							 colData = condition_df,
							 design = ~ condition)
dds

# filtering (default >= 10)
dim(dds)
keep <- rowSums(counts(dds)) >= minCounts
dds <- dds[keep,]
dim(dds)

# designate factor level (control in the front) (optional)
#dds$condition <- factor(dds$condition, levels = c("MRQ","GRD", "DRD"))
# designate control group (optional)
dds$condition <- relevel(dds$condition, ref ="control")
dds$condition

# Use the results() function to analyse differential genes
dds <- DESeq(dds)
# Use the results() function to extract the results of the analysis
#treated vs untreated, log2(treated/untreated)
res <- results(dds)


res <- results(dds, contrast=c("condition","case","control"))
res

# parallel operation
#library("BiocParallel")
#register(MulticoreParam(4))
#dds <- DESeq(dds，parallel = TRUE)

# p-values and adjusted p-values；
# Sort the results in ascending order according to p-value
resOrdered <- res[order(res$pvalue),]
resOrdered

# Number of genes with statistically adjusted p-values < 0.05 (default)
sum(res$padj < pValue, na.rm=TRUE)
# the parameter alpha (adjusted p value cutoff) is 0.1.
#res005 <- results(dds, alpha=0.05)
#summary(res005)
#sum(res05$padj < 0.05, na.rm=TRUE)


# export table of variance analysis results
resOrdered <- resOrdered[complete.cases(resOrdered), ]
resOrdered <- resOrdered[resOrdered$pvalue != 0, ]
write.csv(as.data.frame(resOrdered),
          file=paste0(outFolder,"/diff_genes.csv"), quote = FALSE)


# If only differential genes are exported, the subset function can be used
resSig <- subset(resOrdered, padj < pValue)
write.csv(as.data.frame(resSig),
          file=paste0(outFolder,"/diff_genes_padj",pValue,".csv"),quote = FALSE)

# Visualisation and export of analysis results
# calculate outliers
pdf(paste0(outFolder,"/outliers.pdf"), width=6,height=6)
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
dev.off()

print(paste0("alpha: ",metadata(res)$alpha))
print(metadata(res)$filterThreshold)

if (nrow(condition_df) <= 50) {
  #vsd <- vst(dds, blind=FALSE)
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  rld <- rlog(dds, blind=FALSE)
  
  #log2(n + 1)
  ntd <- normTransform(dds)
  
  head(assay(vsd), 3)
  head(assay(rld), 3)
  head(assay(ntd), 3)
  
  # export normalised data (which can be used for downstream analysis such as GSEA)
  write.csv(as.data.frame(assay(vsd)),
  		  file=paste0(outFolder,"/count_transformation_vst.csv"))
  write.csv(as.data.frame(assay(rld)),
  		  file=paste0(outFolder,"/count_transformation_rlog.csv"))
  write.csv(as.data.frame(assay(ntd)),
  		  file=paste0(outFolder,"/count_transformation_normTransform.csv"))
  
  # distribution visualisation of standardised counts
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
  
# plot heat map of sample distance clustering
# calculate sample distance matrix
sampleDists <- dist(t(assay(vsd)))

# from vector to matrix
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$condition
colnames(sampleDistMatrix) <- vsd$sample

write.csv(as.data.frame(sampleDistMatrix),
          file=paste0(outFolder,"/sample_clusters.csv"))

