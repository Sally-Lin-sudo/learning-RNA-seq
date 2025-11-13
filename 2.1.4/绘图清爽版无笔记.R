# ========== 设置路径 ==========
data_root <- "/Users/Zhuanz/RNA-seq/tumor-transcriptome-demo"
outdir <- file.path(data_root, "2.1.4_output")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ========== 读取文件列表 ==========
files <- c(
  list.files(file.path(data_root, "COAD"), pattern="\\.txt$", full.names=TRUE),
  list.files(file.path(data_root, "READ"), pattern="\\.txt$", full.names=TRUE),
  list.files(file.path(data_root, "ESCA"), pattern="\\.txt$", full.names=TRUE)
)
head(files)

# ========== 合并 counts ==========
counts_list <- list()
for (f in files) {
  df <- read.table(f, header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
  df <- df[, c(1, 7)]                           # geneid 和 count（featureCounts 格式）
  colnames(df) <- c("Geneid", basename(f))
  counts_list[[f]] <- df
}
counts_matrix <- Reduce(function(x, y) merge(x, y, by = "Geneid"), counts_list)
rownames(counts_matrix) <- counts_matrix$Geneid
counts_matrix <- counts_matrix[, -1]            # 删除 Geneid 列
counts_matrix <- as.matrix(counts_matrix)
mode(counts_matrix) <- "numeric"

# ========== 计算 CPM 和 logCPM ==========
CPM_matrix <- t(1e6 * t(counts_matrix) / colSums(counts_matrix))
logCPM_matrix <- log10(CPM_matrix + 1)          # +1 避免 log(0)

# ========== 检查异常并清洗 ==========
# 检查非有限值
cat("non-finite in logCPM:", sum(!is.finite(logCPM_matrix)), "\n")
# 检查 NA
cat("rows with any NA:", sum(rowSums(is.na(logCPM_matrix)) > 0), "\n")
# 检查 sd==0
num_sd0 <- sum(apply(logCPM_matrix, 1, sd) == 0)
cat("rows with sd==0:", num_sd0, "\n")

# 选择“良好”的行（无 NA 且 sd > 0）
sd_row <- apply(logCPM_matrix, 1, sd)
good_rows <- which(sd_row > 0 & rowSums(is.na(logCPM_matrix)) == 0)
length(good_rows)   # 查看保留了多少基因

logCPM_clean <- logCPM_matrix[good_rows, , drop = FALSE]

# ========== 计算 z-score（用清理后的矩阵）==========
z_matrix <- (logCPM_clean - rowMeans(logCPM_clean)) / apply(logCPM_clean, 1, sd)

# 再次确保没有 NA/NaN/Inf
cat("z_matrix non-finite:", sum(!is.finite(z_matrix)), "\n")

# 限制取值范围（clip）
z_matrix[z_matrix > 2] <- 2
z_matrix[z_matrix < -2] <- -2

# ========== 绘图并保存 ==========
library(pheatmap)
# 在屏幕上显示
pheatmap(
  z_matrix,
  show_rownames = FALSE,
  show_colnames = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Tumor transcriptome: logCPM Z-score heatmap"
)

# 保存为 pdf（注意使用 outdir）
pdf(file.path(outdir, "tumor_heatmap.pdf"), width = 8, height = 6)
pheatmap(
  z_matrix,
  show_rownames = FALSE,
  show_colnames = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Tumor transcriptome: logCPM Z-score heatmap"
)
dev.off()
