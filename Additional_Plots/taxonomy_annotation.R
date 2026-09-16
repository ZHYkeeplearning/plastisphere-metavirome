# 加载必要的库
library(readxl)
library(dplyr)
library(ggplot2)
library(ggsci)

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

# 从命令行参数获取输入和输出路径
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_dir <- args[2]
type <- args[3]

# input_file <- "E:/XuYunProject/R绘图项目/物种门水平柱状图/input/gene_catalog_taxonomy_abundance.xls"
# output_dir <- "E:/XuYunProject/R绘图项目/物种门水平柱状图/output/"
# type <- "gene"

print("命令行参数已加载。")

# 读取输入数据
print(paste("正在从文件读取数据：", input_file))
data <- read.csv(input_file, sep = '\t', header = TRUE)

# 显示数据的前几行
print("数据已成功加载。数据前几行预览：")
print(head(data))

# 统计每个物种的条目数量并写入文件
print("正在统计每个物种的条目数量...")
species_stat <- data %>%
  group_by(Species) %>%
  summarise(Count = n(), .groups = "drop") %>%
  arrange(desc(Count))

print("物种统计信息已生成。正在保存到文件...")
write.table(species_stat, paste0(output_dir, "species_stat.xls"), sep = '\t', quote = FALSE, row.names = FALSE)

# 计算每个物种的基因或 contig 数量
print("正在计算每个物种的基因数量...")
species_counts <- data %>%
  group_by(Phylum, Species) %>%
  summarise(Count = n(), .groups = "drop")

# 选择前 8 个门
print("正在根据物种数量选择前 8 个门...")
top_phyla <- species_counts %>%
  group_by(Phylum) %>%
  summarise(total_species = n(), .groups = "drop") %>%
  arrange(desc(total_species)) %>%  # 按 total_species 降序排序
  head(n = 8) %>%                   # 严格筛选出前 8 行
  select(Phylum)

# 选择每个门中基因数量最多的前 5 个物种
print("正在选择每个门中基因数量最多的前 5 个物种...")
top_species_per_phylum <- species_counts %>%
  semi_join(top_phyla, by = "Phylum") %>%
  group_by(Phylum) %>%
  top_n(5, wt = Count) %>%
  arrange(Phylum, desc(Count)) %>%
  ungroup() %>%
  mutate(Species = factor(Species, levels = unique(Species)))

# 过滤掉标记为“unclassified”的物种
# print("正在过滤掉标记为‘unclassified’的物种...")
# top_species_per_phylum <- top_species_per_phylum %>%
#   filter(!grepl("unclassified", Species, ignore.case = TRUE))

print("数据处理完成。正在保存结果到文件...")

out_file = ''
bar_y_title = ''
if(type == 'contig'){
  out_file = "contig_catalog_species_stats.xls"
  bar_y_title = "Number of Contigs"
}else if(type == 'gene'){
  out_file = "gene_catalog_species_stats.xls"
  bar_y_title = "Number of Genes"
}
write.table(top_species_per_phylum, paste0(output_dir, out_file), sep = '\t', quote = FALSE, row.names = FALSE)

# 创建柱状图
print("正在创建柱状图...")
classify_plot <- ggplot(top_species_per_phylum, aes(x = Species, y = Count, fill = Phylum)) +
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
  labs(x = "", y = bar_y_title, title = "Phylum Classification") +
  scale_fill_npg()

print("柱状图已创建。正在显示图表...")
print(classify_plot)

# 将柱状图保存为 PDF 文件
print("正在将柱状图保存为 PDF 文件...")
ggsave(paste0(output_dir, "NR_Species_count.pdf"), plot = classify_plot, width = 11, height = 8, units = "in")

print("分析完成！所有文件已成功保存。")