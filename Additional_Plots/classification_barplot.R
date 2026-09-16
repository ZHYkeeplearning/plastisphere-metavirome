library(tidyverse)
library(ggplot2)
library(scales)
library(ggsci)  # 添加这行以使用 scale_fill_npg()
library(extrafont)
library(showtext)


# 检测操作系统
os <- Sys.info()["sysname"]

# 加载 Times New Roman 字体
if (os == "Windows") {
  # Windows 系统：使用 showtext 加载字体
  font_path <- "C:/Windows/Fonts/times.ttf"  # Times New Roman 字体路径
  if (file.exists(font_path)) {
    font_add("Times New Roman", font_path)
  } else {
    warning("Times New Roman 字体未找到，使用默认字体 Arial。")
    font_add("Arial", "arial.ttf")  # 使用默认字体
  }
} else if (os == "Linux") {
  # Linux 系统：使用 showtext 加载字体
  font_path <- "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf"  # Times New Roman 字体路径
  if (file.exists(font_path)) {
    font_add("Times New Roman", font_path)
  } else {
    warning("Times New Roman 字体未找到，使用默认字体 DejaVu Sans。")
    font_add("DejaVu Sans", "DejaVuSans.ttf")  # 使用默认字体
  }
} else {
  warning("不支持的操作系统，使用默认字体 Arial。")
  font_add("Arial", "arial.ttf")  # 使用默认字体
}

# 启用 showtext
showtext_auto()

# 设置默认字体
if ("Times New Roman" %in% font_families()) {
  main_font <- "Times New Roman"
} else if ("Arial" %in% font_families()) {
  main_font <- "Arial"
} else {
  main_font <- "DejaVu Sans"
}

# 打印使用的字体
cat("使用的字体：", main_font, "\n")


# 设置默认字体
if ("Times New Roman" %in% fonts() || "Times New Roman" %in% font_families()) {
  main_font <- "Times New Roman"
} else {
  main_font <- "DejaVu Sans"
}

# input_file <- "E:/XuYunProject/R绘图项目/分类柱状图/input1/go_categorie.xls"
# output_dir <- "E:/XuYunProject/R绘图项目/分类柱状图/output1/"
# plot_name <- "go"


args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_dir <- args[2]
plot_name <- args[3]

# 读取数据文件
data <- read.csv(input_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# 动态获取列名
col_names <- colnames(data)
level1_col <- col_names[1]
level2_col <- col_names[2]
sample_names <- col_names[3:length(col_names)]

# 获取 Level2 的唯一值，保持原始顺序
level2_order <- unique(data[[level2_col]])

# 创建总和数据
data_sum <- data %>%
  mutate(Sum = rowSums(select(., all_of(sample_names)))) %>%
  select(all_of(c(level1_col, level2_col, "Sum"))) %>%
  mutate(!!level2_col := factor(!!sym(level2_col), levels = level2_order))

# 创建总和图表
p_sum <- ggplot(data_sum, aes(x = !!sym(level2_col), y = Sum, fill = !!sym(level1_col))) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_y_continuous(labels = scales::comma) +
  theme(
    text = element_text(family = main_font),  # 设置全局字体
    axis.text.x = element_text(angle = 65, hjust = 1),
    legend.title = element_blank(),
    legend.position = "right",
    plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "cm"),
    plot.title = element_text(size = 14),  # 主标题
    axis.title = element_text(size = 12),   # 坐标轴标题
    axis.text = element_text(size = 10),    # 坐标轴刻度文本
    legend.text = element_text(size = 10)   # 图例文本
  ) +
  labs(x = "", y = "Number of Genes") +
  scale_fill_npg()

# 保存总和图表
ggsave(paste0(output_dir, "gene_catalog_", plot_name, "_categorie.pdf"), 
       p_sum, width = 12, height = 8, units = "in")

# 为每个样本创建一个图表
for (sample in sample_names) {
  # 准备当前样本的数据
  data_sample <- data %>%
    select(all_of(c(level1_col, level2_col, sample))) %>%
    rename(Count = !!sample) %>%
    mutate(!!level2_col := factor(!!sym(level2_col), levels = level2_order))
  
  # 创建图表
  p <- ggplot(data_sample, aes(x = !!sym(level2_col), y = Count, fill = !!sym(level1_col))) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    scale_y_continuous(labels = scales::comma) +
    theme(
      text = element_text(family = main_font),  # 设置全局字体
      axis.text.x = element_text(angle = 65, hjust = 1),
      legend.title = element_blank(),
      legend.position = "right",
      plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "cm"),
      plot.title = element_text(size = 14),  # 主标题
      axis.title = element_text(size = 12),   # 坐标轴标题
      axis.text = element_text(size = 10),    # 坐标轴刻度文本
      legend.text = element_text(size = 10)   # 图例文本
    ) +
    labs(x = "", y = "Number of Genes") +
    scale_fill_npg()
  
  # 打印图表
  print(p)
  
  # 保存图表
  ggsave(paste0(output_dir, sample, ".gene_", plot_name, "_categorie.pdf"), 
         p, width = 12, height = 8, units = "in", family = main_font)
}
