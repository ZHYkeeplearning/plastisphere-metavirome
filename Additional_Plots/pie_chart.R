library(ggplot2)
library(reshape2)
library(ggsci)
library(dplyr)
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


args <- commandArgs(trailingOnly = TRUE)
legend_title <- args[1]
file_path <- args[2]
output_file <- args[3]

# 
# legend_title="Kingdom"
# file_path <- "E:/XuYunProject/R绘图项目/饼图/多样本_input/Kingdom_abundance.xls"  # 请将此路径替换为你的实际文件路径
# output_file <- "E:/XuYunProject/R绘图项目/饼图/多样本_output/Kingdom_pieplot.pdf"

# 读取数据
data <- read.table(file_path, 
                   sep="\t", 
                   header=TRUE, 
                   check.names=FALSE, 
                   quote="", 
                   stringsAsFactors=FALSE,
                   fileEncoding="UTF-8")

first_col_name <- colnames(data)[1]
numeric_cols <- colnames(data)[-1]
data[numeric_cols] <- lapply(data[numeric_cols], as.numeric)

# 计算每行的总和
data$sum <- rowSums(data[numeric_cols])

# 获取前10行数据
top_data <- data %>%
  arrange(desc(sum)) %>%
  head(10) %>%
  select(-sum)

# 计算其他行的总和
other_sums <- data %>%
  filter(!get(first_col_name) %in% top_data[[first_col_name]]) %>%
  select(all_of(numeric_cols)) %>%
  colSums()

# 创建"Other"行
other_row <- data.frame(stringsAsFactors = FALSE)
other_row[1, first_col_name] <- "Other"
other_row[1, numeric_cols] <- other_sums

# 合并数据
final_data <- rbind(top_data, other_row)

# 准备绘图数据
final_data$Class <- final_data[[first_col_name]]
final_data <- final_data %>% select(-all_of(first_col_name))
melted_data <- melt(final_data, id.vars="Class", variable.name="Sample", value.name="Value")

# 计算样本数量
sample_count <- length(unique(melted_data$Sample))

# 计算合适的行列数
n_col <- ceiling(sqrt(sample_count))
n_row <- ceiling(sample_count / n_col)

# 计算合适的图片尺寸
# 每个饼图的基础大小
base_size <- 2
# 考虑图例和标题的额外空间
extra_width <- 2
extra_height <- 1

width <- n_col * base_size + extra_width
height <- n_row * base_size + extra_height

# 定义颜色
color_palette <- c(pal_npg("nrc")(10), "grey55")

# 绘制饼图
p <- ggplot(melted_data, aes(x="", y=Value, fill=Class)) +
  geom_bar(width=1, stat="identity") +
  coord_polar(theta="y") +
  facet_wrap(~ Sample, scales="free", ncol=n_col) +
  theme_minimal() +
  theme(
    text = element_text(family = main_font),  # 设置全局字体
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    panel.grid=element_blank(),
    strip.text=element_text(size=10),
    plot.title = element_text(size = 14),  # 主标题
    axis.title = element_text(size = 12),   # 坐标轴标题
    legend.title = element_text(size = 10), # 图例标题
    legend.text = element_text(size = 10)   # 图例文本
  ) +
  labs(title="", x="", y="", fill=legend_title) +
  scale_fill_manual(values=color_palette)

# 保存图片
ggsave(output_file, plot=p, width=width, height=height, dpi=300)
