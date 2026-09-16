# 加载必要的库
library(VennDiagram)
library(RColorBrewer)
library(dplyr)
library(extrafont)
loadfonts(device = "pdf")

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_pdf <- args[2]

# input_file <- "E:/XuYunProject/R绘图项目/韦恩图/input/Class_abundance1.xls"
# output_pdf <- "E:/XuYunProject/R绘图项目/韦恩图/output/Class_vennplot.pdf"
out_xls <- gsub("\\.pdf$", ".xls", output_pdf)  # 自动生成Excel文件名

# 读取数据
data_abundance <- read.csv(input_file, sep="\t", header=TRUE, check.names = FALSE, row.names = 1)
print(data_abundance)

# 初始化集合列表
sets <- list()

# 遍历每个样本列，生成不为0的项的集合
for (sample in colnames(data_abundance)) {
  sets[[sample]] <- rownames(data_abundance)[data_abundance[[sample]] != 0]
}

# 准备颜色调色板
num_samples <- length(sets)
if (num_samples == 2) {
  myCol <- c("#E41A1C", "#377EB8")
} else if (num_samples >= 3 && num_samples <= 5) {
  myCol <- brewer.pal(num_samples, "Set1")
} else {
  stop("样本数量超出支持范围（2-5）")
}

# 创建PDF文件
pdf(file = output_pdf, height = 8, width = 8, family = "Times New Roman")

# 创建Venn图 (关键修改部分)
venn.plot <- venn.diagram(
  x = sets,
  category.names = colnames(data_abundance),
  filename = NULL,
  output = TRUE,
  
  # 核心修改参数：固定圆的比例和大小
  scaled = FALSE,        # 禁用自动比例缩放
  radius = rep(1, num_samples),  # 所有圆相同半径
  
  # 圆圈样式
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # 数字样式
  cex = 1.2,  # 适当放大数字
  fontfamily = "Times New Roman",  # 设置字体
  
  # 集合名称样式
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.fontfamily = "Times New Roman",  # 分类名称字体
  
  # 其他布局调整
  margin = 0.15  # 增加边距防止文字截断
)

# 绘制Venn图
grid.draw(venn.plot)

# 关闭PDF设备
dev.off()



# 计算最大集合长度
max_length <- max(sapply(sets, length))

# 用空字符串填充每个集合，使得它们长度一致
sets_padded <- lapply(sets, function(x) {
  length(x) <- max_length
  # 如果长度不足，填充为空字符串
  x[is.na(x)] <- ""
  return(x)
})

# 将列表转换为数据框
sets_df <- as.data.frame(sets_padded)

# 使用 write.table 保存为 xls 文件
write.table(sets_df, file = out_xls, sep = "\t", quote = FALSE, row.names = FALSE)
