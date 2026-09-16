# --- 步骤 0: 安装并加载必要的包 ---
# 确保已安装这些包
# install.packages(c("tidyverse", "vegan", "ggpubr", "rstatix", "readr", "rcompanion", "showtext"))
library(tidyverse)
library(vegan)
library(ggpubr)
library(rstatix)
library(readr)
library(rcompanion)
library(showtext)

# --- 步骤 1, 2, 3, 4: 数据读取与处理 (与之前相同) ---
abundance_file <- "input/contig_count_abundance.xls"
metadata_file  <- "input/sample_group.tsv"
output_dir     <- "output"

cat("--> 读取数据...\n")
abundance_data <- read.delim(abundance_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
metadata <- read.delim(metadata_file, sep = "\t", header = TRUE, row.names = 1)
if ("Detail" %in% colnames(abundance_data)) {
  abundance_data <- abundance_data %>% dplyr::select(-Detail)
}
community_data <- t(abundance_data)
community_data <- community_data[rownames(metadata), , drop = FALSE]
community_data_num <- matrix(as.numeric(community_data), nrow = nrow(community_data), ncol = ncol(community_data))
rownames(community_data_num) <- rownames(community_data)
colnames(community_data_num) <- colnames(community_data)

cat("--> 计算Alpha多样性...\n")
alpha_chao <- apply(community_data_num, 1, function(x) tryCatch(vegan::estimateR(x)["S.chao1"], error = function(e) NA))
alpha_shannon <- vegan::diversity(community_data_num, index = "shannon")
alpha_simpson <- vegan::diversity(community_data_num, index = "simpson")
alpha_invsimpson <- 1 / alpha_simpson
alpha_invsimpson[is.infinite(alpha_invsimpson)] <- NA
alpha_diversity_df <- data.frame(Chao1 = alpha_chao, Shannon = alpha_shannon, Simpson = alpha_simpson, InvSimpson = alpha_invsimpson, row.names = rownames(community_data_num))

cat("--> 合并元数据...\n")
metadata_prepared <- metadata %>% rownames_to_column(var = "Sample")
colnames(metadata_prepared)[2] <- "SampleGroup"
alpha_with_metadata_wide <- alpha_diversity_df %>% rownames_to_column(var = "Sample") %>% left_join(metadata_prepared, by = "Sample")
alpha_long_format <- alpha_with_metadata_wide %>% pivot_longer(cols = c(Chao1, Shannon, Simpson, InvSimpson), names_to = "Diversity_Index", values_to = "Value")



# -------------------------------------------------------------------------- #
# --- 步骤 5 (Kruskal+Dunn+cldList 版): 强制按均值排序字母 ---
# -------------------------------------------------------------------------- #
# 确保加载了 rcompanion
# install.packages("rcompanion")
library(rcompanion)

cat("--> 进行 Kruskal-Wallis + Dunn 检验并生成字母 (已优化排序)...\n")

stat.letters <- data.frame()
all_pairwise_p_values <- data.frame() 

for (index in unique(alpha_long_format$Diversity_Index)) {
  
  # 1. 提取并清洗数据
  subset_data <- dplyr::filter(alpha_long_format, Diversity_Index == index)
  
  # 筛选有效分组（至少2个样本）
  valid_groups <- subset_data %>% 
    dplyr::filter(!is.na(Value)) %>% 
    group_by(SampleGroup) %>% 
    summarise(n = n(), .groups = 'drop') %>% 
    dplyr::filter(n >= 2) %>% 
    pull(SampleGroup)
  
  if (length(valid_groups) < 2) {
    next
  }
  
  subset_data_for_test <- dplyr::filter(subset_data, SampleGroup %in% valid_groups)
  
  # =========================================================================
  # 【关键修改】: 在检验前，将分组因子按均值【降序】排列
  # 这样 dunn_test 的比较顺序就是 "最大值组 vs 最小值组"
  # cldList 就会把 "a" 赋予排列在最前面的组（即最大值组）
  # =========================================================================
  subset_data_for_test$SampleGroup <- factor(
    subset_data_for_test$SampleGroup, 
    levels = levels(reorder(subset_data_for_test$SampleGroup, subset_data_for_test$Value, mean, decreasing = TRUE))
  )
  
  # 2. Kruskal-Wallis 检验
  kruskal_test_res <- kruskal.test(Value ~ SampleGroup, data = subset_data_for_test)
  
  if (kruskal_test_res$p.value < 0.05) {
    # 3. Dunn's Test (Post-hoc)
    dunn_test_res <- dunn_test(data = subset_data_for_test, Value ~ SampleGroup, p.adjust.method = "fdr")
    
    # 保存P值
    p_values_to_save <- dunn_test_res %>% mutate(Diversity_Index = index)
    all_pairwise_p_values <- rbind(all_pairwise_p_values, p_values_to_save)
    
    # 准备 cldList 需要的输入
    p_values <- dunn_test_res$p.adj
    comparison_names <- paste(dunn_test_res$group1, dunn_test_res$group2, sep = "-")
    names(p_values) <- comparison_names
    
    # 4. 生成字母 (cldList)
    # 因为前面重排序了因子，这里生成的字母通常就是有序的 (a=最大)
    letters <- rcompanion::cldList(p.value = p_values, comparison = comparison_names, threshold = 0.05)
    
    # 整理结果
    letters <- letters %>% 
      rename(SampleGroup = Group) %>% 
      mutate(Diversity_Index = index) %>% 
      dplyr::select(SampleGroup, Letter, Diversity_Index)
    
    stat.letters <- rbind(stat.letters, letters)
    
  } else {
    # 如果整体差异不显著，所有组标记为 "a"
    all_groups <- levels(subset_data_for_test$SampleGroup) # 注意用 levels 获取排序后的组
    letters <- data.frame(SampleGroup = all_groups, Letter = "a", Diversity_Index = index)
    stat.letters <- rbind(stat.letters, letters)
  }
}

# 计算字母在图上的Y轴位置 (最大值上方10%)
letter_positions <- alpha_long_format %>% 
  group_by(Diversity_Index, SampleGroup) %>% 
  summarise(y.position = max(Value, na.rm = TRUE) * 1.15, .groups = 'drop')

# 合并
plot_labels <- left_join(stat.letters, letter_positions, by = c("Diversity_Index", "SampleGroup"))

cat("--> 显著性字母生成完毕，准备绘图...\n")

# --- 步骤 6: 准备绘图的通用设置 (字体，颜色等) ---
font_add(family = "TNR", regular = "times.ttf", bold = "timesbd.ttf", italic = "timesi.ttf", bolditalic = "timesbi.ttf")
showtext_auto()
all_groups <- levels(factor(alpha_long_format$SampleGroup))
plot_colors <- c("#74B3C9", "#E68674", "#A6D554", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462")
if (length(all_groups) > length(plot_colors)) { plot_colors <- rainbow(length(all_groups)) }
named_colors <- setNames(plot_colors[1:length(all_groups)], all_groups)

# --- 步骤 7a: 绘制并保存组合图 ---
cat("--> 正在生成组合图...\n")
combined_plot <- ggplot(alpha_long_format, aes(x = SampleGroup, y = Value, fill = SampleGroup)) +
  geom_boxplot(width = 0.6, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  facet_wrap(~ Diversity_Index, scales = "free_y", ncol = 2) +
  geom_text(data = plot_labels, aes(x = SampleGroup, y = y.position, label = Letter), family = "TNR", size = 5, fontface = "bold", color = "black") +
  scale_fill_manual(values = named_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  theme_classic(base_size = 14, base_family = "TNR") +
  labs(title = "Alpha Diversity Comparison", x = "Group", y = "Diversity Index Value") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), plot.subtitle = element_text(hjust = 0.5, size = 10), axis.text = element_text(color = "black", face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), strip.text = element_text(size = 14, face = "bold"), strip.background = element_rect(fill = "grey90", color = NA), legend.position = "none", panel.spacing = unit(1.5, "lines"))

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
plot_width_combined <- 12
plot_height_combined <- 10
output_filename_png <- file.path(output_dir, "Alpha_Diversity_Combined.png")
ggsave(filename = output_filename_png, plot = combined_plot, width = plot_width_combined, height = plot_height_combined, dpi = 300)
output_filename_pdf <- file.path(output_dir, "Alpha_Diversity_Combined.pdf")
ggsave(filename = output_filename_pdf, plot = combined_plot, width = plot_width_combined, height = plot_height_combined)
print(combined_plot)
cat(paste0("\n组合图已保存至: ", output_filename_png, " 和 ", output_filename_pdf, "\n"))

# --- 【新增功能】步骤 7b: 循环绘制并保存每个指数的独立PDF图 ---
cat("--> 正在为每个指数生成独立的PDF文件...\n")
indices_to_plot <- unique(alpha_long_format$Diversity_Index)
plot_width_single <- max(8, length(all_groups) * 0.8) # 为单图计算合适的宽度
plot_height_single <- 7

for (current_index in indices_to_plot) {
  
  # 筛选当前指数的数据和标签
  data_single <- dplyr::filter(alpha_long_format, Diversity_Index == current_index)
  labels_single <- dplyr::filter(plot_labels, Diversity_Index == current_index)
  
  # 创建独立的图表对象 (无 facet_wrap)
  single_plot <- ggplot(data_single, aes(x = SampleGroup, y = Value, fill = SampleGroup)) +
    geom_boxplot(width = 0.6, alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
    geom_text(data = labels_single, aes(x = SampleGroup, y = y.position, label = Letter), family = "TNR", size = 6, fontface = "bold", color = "black") + # 稍微增大了字母
    scale_fill_manual(values = named_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    theme_classic(base_size = 16, base_family = "TNR") + # 稍微增大了基础字号
    labs(
      title = paste("Alpha Diversity -", current_index), # 动态生成标题
      subtitle = "Groups sharing a letter are not significantly different (Dunn's test, p.adj < 0.05)",
      x = "Group",
      y = "Diversity Index Value"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.text = element_text(color = "black", face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "none"
    )
  
  # 构建文件名并保存
  output_filename_single <- file.path(output_dir, paste0("Alpha_Diversity_", current_index, ".pdf"))
  ggsave(filename = output_filename_single, plot = single_plot, width = plot_width_single, height = plot_height_single)
  
  cat(paste("  - 已保存:", output_filename_single, "\n"))
}

# --- 步骤 8: 保存所有数据文件 (与之前相同) ---
showtext_auto(FALSE)
cat("\n--> 正在保存数据文件...\n")
csv_data_filename <- file.path(output_dir, "alpha_diversity_with_metadata.csv")
write.csv(alpha_with_metadata_wide, csv_data_filename, row.names = FALSE, na = "")
cat(paste("多样性指数宽格式表格已保存至:", csv_data_filename, "\n"))
csv_stat_filename <- file.path(output_dir, "alpha_diversity_statistics_with_Letters.csv")
write.csv(stat.letters, csv_stat_filename, row.names = FALSE, na = "")
cat(paste("显著性字母统计结果已保存至:", csv_stat_filename, "\n"))
csv_p_values_filename <- file.path(output_dir, "alpha_diversity_pairwise_p_values.csv")
write.csv(all_pairwise_p_values, csv_p_values_filename, row.names = FALSE, na = "")
cat(paste("详细的两两比较p值结果已保存至:", csv_p_values_filename, "\n"))