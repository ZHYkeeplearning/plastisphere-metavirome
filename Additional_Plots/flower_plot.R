library(plotrix)
library(RColorBrewer)
library(extrafont)
loadfonts(device = "pdf")

# 定义 flower_plot 函数
flower_plot <- function(sample, value, start, a, b, r, overlap,
                        ellipse_col = color.scale(value, alpha=0.2, color.spec= "rgb"),
                        circle_col = "green", circle_text_cex = 1) {
  par(
    bty = "n",
    ann = F,
    xaxt = "n",
    yaxt = "n",
    mar = c(1, 1, 1, 1)
  )
  plot(c(0, 10), c(0, 10), type = "n")
  
  n   <- length(sample)
  deg <- 360 / n
  
  # 动态计算样本名称的文字大小
  sample_text_cex <- min(max(1.2 - (n / 40), 0.4), 1.0)
  
  # 动态计算数字的文字大小
  value_text_cex <- min(max(1.5 - (n / 50), 0.5), 0.9)
  
  # 先绘制花瓣
  lapply(1:n, function(t) {
    draw.ellipse(
      x = 5 + cos((start + deg * (t - 1)) * pi / 180),
      y = 5 + sin((start + deg * (t - 1)) * pi / 180),
      col = ellipse_col[t],
      border = ellipse_col[t],
      a = a,
      b = b,
      angle = deg * (t - 1)
    )
  })
  
  # 再绘制文字
  lapply(1:n, function(t) {
    # 绘制数字，使用动态大小
    text(x = 5 + 2.5 * cos((start + deg * (t - 1)) * pi / 180),
         y = 5 + 2.5 * sin((start + deg * (t - 1)) * pi / 180),
         value[t],
         cex = value_text_cex)
    
    if (deg * (t - 1) < 180 && deg * (t - 1) > 0) {
      # 绘制样本名，使用动态大小
      text(
        x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
        y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
        sample[t],
        srt = deg * (t - 1) - start,
        adj = 1,
        cex = sample_text_cex
      )
      
    } else {
      # 绘制样本名，使用动态大小
      text(
        x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
        y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
        sample[t],
        srt = deg * (t - 1) + start,
        adj = 0,
        cex = sample_text_cex
      )
    }
  })
  
  # 绘制中心圆
  draw.circle(
    x = 5,
    y = 5,
    r = r,
    col = circle_col,
    border = circle_col
  )
  text(x = 5, y = 5, paste('Core: ', overlap))
}

# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
out_pdf <- args[2]
# 
# input_file <- "E:/XuYunProject/R绘图项目/韦恩图/input_63/Class_abundance.xls"
# out_pdf <- "E:/XuYunProject/R绘图项目/韦恩图/output_63/Class_flowerplot1.pdf"


# 自动生成Excel文件名
out_xls <- gsub("\\.pdf$", ".xls", out_pdf)

# 读取数据
data_abundance <- read.csv(input_file, sep="\t", header=TRUE, check.names = FALSE, row.names = 1)
print(data_abundance)

# 初始化集合列表
sets <- list()
for (sample in colnames(data_abundance)) {
  sets[[sample]] <- rownames(data_abundance)[data_abundance[[sample]] != 0]
}

# 取sets的key作为一个数组
B <- names(sets)

# 计算每个集合中的元素数量
B.data <- sapply(sets, length)

# 定义颜色
num_colors <- length(B.data)
color_palette <- colorRampPalette(brewer.pal(8, "Set3"))
B.col <- color_palette(num_colors)

# 计算所有集合的交集
shared_elements <- Reduce(intersect, sets)
shared_elements_count <- length(shared_elements)

# 打印交集结果
print(shared_elements_count)

# 计算减去交集后的值
adjusted_values <- B.data - shared_elements_count
adjusted_values <- pmax(adjusted_values, 0)  # 确保值不小于0

# 使用示例：生成PDF
pdf(out_pdf, family = "Times New Roman")
flower_plot(B, adjusted_values, 90, 0.5, 2, 1, shared_elements_count, ellipse_col = B.col)
dev.off()

# 计算最大集合长度
max_length <- max(sapply(sets, length))

# 用空字符串填充每个集合，使得它们长度一致
sets_padded <- lapply(sets, function(x) {
  length(x) <- max_length
  x[is.na(x)] <- ""
  return(x)
})

# 将列表转换为数据框
sets_df <- as.data.frame(sets_padded)

# 使用 write.table 保存为 xls 文件
write.table(sets_df, file = out_xls, sep = "\t", quote = FALSE, row.names = FALSE)