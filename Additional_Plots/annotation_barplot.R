# 加载必要的库
library(readxl)
library(ggplot2)
library(dplyr)
library(ggsci)
library(tidyr)
library(extrafont)
library(showtext)

# 检测操作系统
os <- Sys.info()["sysname"]

# 加载 Times New Roman 字体
if (os == "Windows") {
  font_path <- "C:/Windows/Fonts/times.ttf"
  if (file.exists(font_path)) {
    font_add("Times New Roman", font_path)
  } else {
    warning("Times New Roman 字体未找到，使用默认字体 Arial。")
    font_add("Arial", "arial.ttf")
  }
} else if (os == "Linux") {
  font_path <- "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf"
  if (file.exists(font_path)) {
    font_add("Times New Roman", font_path)
  } else {
    warning("Times New Roman 字体未找到，使用默认字体 DejaVu Sans。")
    font_add("DejaVu Sans", "DejaVuSans.ttf")
  }
} else {
  warning("不支持的操作系统，使用默认字体 Arial。")
  font_add("Arial", "arial.ttf")
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

cat("使用的字体：", main_font, "\n")

# 主程序
args <- commandArgs(trailingOnly = TRUE)
plot_name <- args[1]
input_file <- args[2]
output_dir <- args[3]
top_n <- as.numeric(args[4])
angle <- as.numeric(args[5])
hjust <- as.numeric(args[6])
x_label_show_type <- as.numeric(args[7])

# 读取数据
data <- read.csv(input_file, sep = '\t', header = TRUE, check.names = FALSE)

# 动态获取第一列的名称
class_column <- names(data)[1]

# 添加处理后的x轴显示列
if(x_label_show_type == 1){
  data <- data %>%
    mutate(Display_Name = sub(":.*", "", !!sym(class_column)))
} else if(x_label_show_type == 2){
  data <- data %>%
    mutate(Display_Name = sub("^(.*?)\\s.*$", "\\1", !!sym(class_column)))
}

# 将数据转换为长格式
data_long <- pivot_longer(data, cols = -c(!!sym(class_column), Display_Name), 
                          names_to = "Sample", values_to = "Value")

# 获取样本列表
samples <- unique(data_long$Sample)

# 定义通用主题
my_theme <- theme_minimal() +
  theme(
    text = element_text(family = main_font),
    plot.background = element_rect(fill = "white", colour = "black"),
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.text.x = element_text(angle = angle, hjust = hjust),
    legend.position = "right",
    legend.box = "vertical",
    legend.key.size = unit(0.5, "cm"),
    legend.title = element_text(size = 10),
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10)
  )

# 为每个样本绘制柱状图
for (sample in samples) {
  sample_data <- filter(data_long, Sample == sample) %>%
    ungroup() %>%
    slice_max(order_by = Value, n = top_n, with_ties = FALSE)
  
  categories <- unique(sample_data[[class_column]])
  col_vector <- colorRampPalette(pal_npg("nrc")(9))(length(categories))
  names(col_vector) <- categories
  
  p <- ggplot(sample_data, aes(x = Display_Name, y = Value, fill = !!sym(class_column))) +
    geom_bar(stat = "identity", position = position_dodge(), colour = "black", alpha = 0.8) +
    scale_fill_manual(values = col_vector) +
    scale_y_continuous(limits = c(0, max(sample_data$Value) * 1.1), expand = c(0, 0)) +
    labs(x = class_column, y = "Gene Number") +
    my_theme +
    guides(fill = guide_legend(ncol = 1))
  
  print(p)
  ggsave(paste(output_dir, sample, '.gene_', plot_name, "_categorie.pdf", sep = ""), 
         plot = p, width = 10, height = 6, dpi = 300)
}

# 总和的柱状图处理
total_values <- data_long %>%
  group_by(!!sym(class_column), Display_Name) %>%
  summarise(Total_Value = sum(Value, na.rm = TRUE), .groups = "drop") %>%
  slice_max(order_by = Total_Value, n = top_n, with_ties = FALSE)

filtered_categories <- unique(total_values[[class_column]])
filtered_col_vector <- colorRampPalette(pal_npg("nrc")(9))(length(filtered_categories))
names(filtered_col_vector) <- filtered_categories

p_total <- ggplot(total_values, aes(x = Display_Name, y = Total_Value, fill = !!sym(class_column))) +
  geom_bar(stat = "identity", position = position_dodge(), colour = "black", alpha = 0.8) +
  scale_fill_manual(values = filtered_col_vector) +
  scale_y_continuous(limits = c(0, max(total_values$Total_Value) * 1.1), expand = c(0, 0)) +
  labs(x = class_column, y = "Gene Number") +
  my_theme +
  guides(fill = guide_legend(ncol = 1))

print(p_total)
ggsave(paste0(output_dir, "gene_catalog_", plot_name, "_categorie.pdf"), 
       plot = p_total, width = 10, height = 6, dpi = 300)