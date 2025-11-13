###
# 0.准备工作
# install.packages(c("pheatmap","matrixStats"))
# 如果想用 edgeR 计算 CPM，可用 BiocManager 安装：
# if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
# BiocManager::install("edgeR")

###
# 1. 设置工作环境、设置数据的根目录&结果输出目录
data_root <- "/Users/Zhuanz/RNA-seq/tumor-transcriptome-demo"
outdir <- file.path(data_root, "2.1.4_output")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

###
#2.跳过递归查找，已知全部都是feature count，直接读取文件
files <- c(
  list.files(file.path(data_root, "COAD"), pattern="\\.txt$", full.names=TRUE),
  list.files(file.path(data_root, "READ"), pattern="\\.txt$", full.names=TRUE),
  list.files(file.path(data_root, "ESCA"), pattern="\\.txt$", full.names=TRUE)
)
# 查看前几个文件，确保路径正确 
head(files)

###
# 3. 读取并合并为一个表达矩阵
# 定义一个空列表用于存放每个样本的 counts 
counts_list <- list()
for (f in files) {
  df <- read.table(f, header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
  
  # 提取 gene_id 和 counts（featureCounts 输出第 1 列为 gene_id，第 7 列为 count）
  df <- df[, c(1, 7)]
  colnames(df) <- c("Geneid", basename(f))
  
  # 存入列表
  counts_list[[f]] <- df 
}

# 以 Geneid 为基准合并所有样本 
counts_matrix <- Reduce(function(x, y) merge(x, y, by = "Geneid"), counts_list)

# 将 Geneid 设为行名
rownames(counts_matrix) <- counts_matrix$Geneid
# 去掉第一列 Geneid
counts_matrix <- counts_matrix[, -1] 

# 检查合并结果
dim(counts_matrix) # 行为基因数，列为样本数
head(counts_matrix[, 1:5])

######
# 4. 计算 logCPM（每百万标准化后取 log10）
######

# 将矩阵转换为数值型
counts_matrix <- as.matrix(counts_matrix)
mode(counts_matrix) <- "numeric"

# CPM = 每个基因的 count / 每个样本的总 reads × 1,000,000
CPM_matrix <- t(1e6 * t(counts_matrix) / colSums(counts_matrix))
logCPM_matrix <- log10(CPM_matrix + 1)  # 加上 pseudocount 1，避免 log(0)

## 或者利用R包edgeR的cpm函数计算也是可以的:
## 如果没有安装edgeR，可通过BiocManager::install("edgeR")安装
# library(edgeR)
# count.matrix 行为基因，列为样本，数值为counts
#y <- DGEList(counts = count.matrix) # 定义edgeR用于存储基因表达信息的DGEList对象
#CPM.matrix <- edgeR::cpm(y,log=F) # 计算CPM
#log10.CPM.matrix <- log10(CPM.matrix+1) # 1 为pseudocount, 避免count为0时对数未定义的情况 

# ========== 检查异常并清洗 ==========
# 检查非有限值
cat("non-finite in logCPM:", sum(!is.finite(logCPM_matrix)), "\n")
# 检查 NA
cat("rows with any NA:", sum(rowSums(is.na(logCPM_matrix)) > 0), "\n")
# 检查 sd==0
num_sd0 <- sum(apply(logCPM_matrix, 1, sd) == 0)
cat("rows with sd==0:", num_sd0, "\n")   # 发现有 103个😮

# 选择goodrows（无 NA 且 sd > 0）
sd_row <- apply(logCPM_matrix, 1, sd)
good_rows <- which(sd_row > 0 & rowSums(is.na(logCPM_matrix)) == 0)
length(good_rows)   # 查看保留了多少基因：1897😮

logCPM_clean <- logCPM_matrix[good_rows, , drop = FALSE]


# 计算Z：每个基因（行）在不同样本中的 Z score
z_matrix <- (logCPM_clean - rowMeans(logCPM_clean)) / apply(logCPM_clean, 1, sd)
  # apply(log10.CPM.matrix,1,sd)表示计算每行(1表示行,2表示列)的标准差(sd函数)
  # rowMeans(log10.CPM.matrix)和apply(log10.CPM.matrix,1,mean)效果是一样的

# 再次确保没有 NA/NaN/Inf
cat("z_matrix non-finite:", sum(!is.finite(z_matrix)), "\n")

# clip [-2, 2]，避免颜色过亮 
z_matrix[z_matrix > 2] <- 2 
z_matrix[z_matrix < -2] <- -2

# z_matrix <- na.omit(z_matrix) # 去掉含 NA 的行


# ========== 绘图并保存 ==========
# if (!require(pheatmap)) install.packages("pheatmap") 
library(pheatmap)

pheatmap(
  z_matrix,
  show_rownames = FALSE, # 不显示基因名（太多）
  show_colnames = FALSE, # 不显示样本名
  color = colorRampPalette(c("blue", "white", "red"))(100),
  clustering_distance_cols = "euclidean", # 样本间距离
  clustering_method = "complete",
  main = "Tumor transcriptome: logCPM Z-score heatmap" 
)

# 保存结果为pdf
pdf(file.path(outdir, "tumor_heatmap.pdf"), width = 8, height = 6)
pheatmap(
  z_matrix,
  show_rownames = FALSE, # 不显示基因名（太多）
  show_colnames = FALSE, # 不显示样本名
  color = colorRampPalette(c("blue", "white", "red"))(100),
  clustering_distance_cols = "euclidean", # 样本间距离
  clustering_method = "complete",
  main = "Tumor transcriptome: logCPM Z-score heatmap" 
)
dev.off()

