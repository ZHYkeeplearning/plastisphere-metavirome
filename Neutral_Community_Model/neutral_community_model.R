# ===================================================================
# 步骤一：加载包
# ===================================================================
# 清理工作环境
rm(list = ls())

# 加载所有需要的库
library(MicEco)
library(ggplot2)
library(ggimage)
library(plyr)
input_file <- "input/abundance.tsv"
output_dir <- "output"
output_prefix <- "vOTU"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)



# ===================================================================
# 步骤二：数据读取与稳健的预处理
# ===================================================================

# 1. 使用正确的方式读取TSV文件
read_table <- read.table(file = input_file, 
                         row.names = 1,      # 第一列是行名
                         header = TRUE,      # 文件有表头
                         sep = "\t",         # 分隔符是制表符 (Tab)
                         check.names = FALSE) # 防止R自动修改不规范的列名

# 2. 过滤掉总丰度为0的物种（行）
cat("原始vOTU数量: ", nrow(read_table), "\n")
read_table_filt1 <- read_table[rowSums(read_table) > 0, ]
cat("过滤稀有vOTU后数量: ", nrow(read_table_filt1), "\n")

# 3. 过滤掉总丰度为0的样本（列）
cat("原始样本数量: ", ncol(read_table_filt1), "\n")
read_table_filt2 <- read_table_filt1[, colSums(read_table_filt1) > 0]
cat("过滤空样本后数量: ", ncol(read_table_filt2), "\n")

# 4. 将清理后的数据转置为 "样本 x 物种" 格式
data_for_fit <- t(read_table_filt2)


# ===================================================================
# 步骤三：运行中性模型分析
# ===================================================================

# 使用清理并转置好的数据运行模型
res <- neutral.fit(data_for_fit)

# 提取模型结果
m <- res[[1]][1]
N <- res[[1]][4]
Nm <- N * m
r2 <- res[[1]][3]
out <- res[[2]]

# ===================================================================
# 步骤四：绘图 (已修正饼图垂直位置)
# ===================================================================

# 处理数据以进行绘图 (这部分不变)
out$group <- with(out, ifelse(freq < Lower, "#509579",
                              ifelse(freq > Upper, "#cf9198", "#485970")))

# 绘制模型结果图 (p1的基础部分不变)
p1 <- ggplot(data = out) +
  geom_line(aes(x = log(p), y = freq.pred), linewidth = 1.2, linetype = 1) +
  geom_line(aes(x = log(p), y = Lower), linewidth = 1.2, linetype = 2) +
  geom_line(aes(x = log(p), y = Upper), linewidth = 1.2, linetype = 2) +
  geom_point(aes(x = log(p), y = freq, color = group), size = 2) +
  xlab("log10(mean relative abundance)") +
  ylab("Occurrence frequency") +
  scale_colour_manual(values = c("#485970", "#cf9198", "#509579")) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(linewidth = 0.5, colour = "black"),
    axis.line.y = element_line(linewidth = 0.5, colour = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 18),
    legend.position = "none",
    text = element_text(family = "sans", size = 18),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# 计算并准备饼图数据 (这部分不变)
data_summary <- as.data.frame(table(out$group))
colnames(data_summary) <- c("group", "nums")
data_summary$group <- factor(data_summary$group, levels = c("#485970", "#509579", "#cf9198"))
data_summary <- data_summary[order(data_summary$group), ]
data_summary$type <- c("med", "low", "high")
data_summary$percentage <- round(data_summary$nums / sum(data_summary$nums) * 100, 1)
data_summary$label <- paste(data_summary$type, paste(data_summary$percentage, "%", sep = ''))

p2 <- ggplot(data_summary, aes(x = "", y = nums, fill = group)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(
    values = c("#485970", "#509579", "#cf9198"),
    labels = data_summary$label
  ) +
  theme_void() +
  theme(legend.text = element_text(size = 12))


# --- 定位注释元素 ---
# 1. R2, Nm, m 文本 (使用Inf/-Inf精确定位到角落，这部分不变)
p1_with_text <- p1 +
  annotate("text", x = Inf, y = -Inf, label = paste("R2 =", round(r2, 3)), size = 5, hjust = 1.1, vjust = -3.5) +
  annotate("text", x = Inf, y = -Inf, label = paste("Nm =", round(Nm, 0)), size = 5, hjust = 1.1, vjust = -2.0) 
# +
#   annotate("text", x = Inf, y = -Inf, label = paste("m =", round(m, 3)), size = 5, hjust = 1.1, vjust = -0.5)

# --- 关键修改：在这里重新计算饼图的位置 ---
# 2.A 获取X轴范围，用于计算饼图的宽度和水平位置 (这部分不变)
x_range <- range(log(out$p), na.rm = TRUE)
pie_width <- (x_range[2] - x_range[1]) * 0.3 # 宽度占总宽度的30%
pie_xmin <- x_range[1]
pie_xmax <- x_range[1] + pie_width

# 2.B 手动指定饼图的高度和垂直位置 (关键改动)
#     您可以微调下面这两个值，来控制饼图的大小和位置
pie_height_ratio <- 0.35  # 将饼图的高度设置为Y轴可见范围(0-1)的35%
pie_top_y        <- 1.05  # 将饼图的顶部放在Y=1.05的位置 (略高于1.0，确保贴顶)

# 根据上面的参数计算最终的Y坐标
pie_height <- (1 - 0) * pie_height_ratio # Y轴可见范围是1，所以高度就是比例值
pie_ymax <- pie_top_y
pie_ymin <- pie_top_y - pie_height


# 3. 组合所有元素 (这部分不变)
p_final <- p1_with_text +
  annotation_custom(
    grob = ggplotGrob(p2),
    xmin = pie_xmin,
    xmax = pie_xmax,
    ymin = pie_ymin,
    ymax = pie_ymax
  ) +
  coord_cartesian(clip = 'off') +
  theme(plot.margin = ggplot2::margin(t = 1, r = 1, b = 0.5, l = 0.5, unit = "cm"))

# ===================================================================
# 步骤五：显示并保存图形
# ===================================================================

# 在RStudio中显示最终图形
print(p_final)

# --- 新增代码：将图形保存为PDF文件 ---
ggsave(
  file.path(output_dir, paste0(output_prefix, "_Neutral_Community_Model_Output.pdf")), # 1. 文件名：你可以自定义文件名
  plot = p_final,                     # 2. 要保存的图形对象
  width = 10,                         # 3. 图像的宽度（单位：英寸）
  height = 8                          # 4. 图像的高度（单位：英寸）
)

cat("图形已成功保存为 Neutral_Community_Model_Output.pdf\n")

# ===================================================================
# 步骤六：保存结果表格
# ===================================================================

# --- 新增代码：将关键的结果数据框保存为文件 ---

# 1. 保存详细的vOTU拟合结果 (out 表格)
write.table(
  out,
  file = file.path(output_dir, paste0(output_prefix, "_Neutral_Model_Detailed_Results.tsv")),
  sep = "\t",          # 使用制表符分隔，保存为tsv文件
  quote = FALSE,       # 不给字符串加引号，方便后续程序读取
  row.names = TRUE     # 保留行名，因为行名是vOTU的ID
)

# 2. 保存汇总统计结果 (data_summary 表格)
write.table(
  data_summary,
  file = file.path(output_dir, paste0(output_prefix, "_Neutral_Model_Summary_Stats.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE    # 汇总表的行号没有特殊意义，无需保留
)

cat("必要的表格已成功保存到输出文件夹中。\n")