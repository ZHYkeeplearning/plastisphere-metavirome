library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggsci)
library(tidyr)
library(alluvial)
library(ggalluvial)
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

plot_abundance <- function(type, input_file, group_file, top_n, output_folder) {
  
  ################################ 数据预处理 ##################################
  
  # 读取并处理数据
  data_sample <- tryCatch({
    # 读取数据
    raw_data <- read.csv(input_file, sep = '\t', header = TRUE, 
                         check.names = FALSE, 
                         stringsAsFactors = FALSE)
    
    # 移除全为NA的行和列
    raw_data <- raw_data[rowSums(!is.na(raw_data)) > 0, ]
    raw_data <- raw_data[, colSums(!is.na(raw_data)) > 0]
    
    # 移除值为0或NA的行
    # 确保 raw_data 是数据框
    if (!is.data.frame(raw_data)) {
      raw_data <- as.data.frame(raw_data)
    }
    
    # 检查 numeric_cols
    numeric_cols <- sapply(raw_data, is.numeric)
    
    # 确保至少有一个数值列
    if (sum(numeric_cols) == 0) {
      stop("数据中没有数值列，请检查数据读取步骤。")
    }
    
    # 计算行和
    row_sums <- rowSums(raw_data[, numeric_cols, drop = FALSE], na.rm = TRUE)
    
    raw_data <- raw_data[row_sums > 0, ]
    
    first_col_name <- names(raw_data)[1]
    
    # 将除第一列外的所有列转换为数值型，并将NA转换为0
    for(col in names(raw_data)[-1]) {
      raw_data[[col]] <- as.numeric(as.character(raw_data[[col]]))
      raw_data[[col]][is.na(raw_data[[col]])] <- 0
    }
    
    raw_data
  }, error = function(e) {
    stop("读取输入文件时出错: ", e$message)
  })
  
  first_col_name <- names(data_sample)[1]
  
  # 转换为长格式数据
  data <- data_sample %>%
    pivot_longer(
      cols = -first_col_name, 
      names_to = "Sample", 
      values_to = "Value"
    ) %>%
    select(Sample, Category = first_col_name, Value) %>%
    # 移除NA值和0值
    filter(!is.na(Value), Value > 0)
  
  # 计算排名前N的分类
  category_sums <- data %>%
    group_by(Category) %>%
    summarise(TotalValue = sum(Value, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(TotalValue)) %>%
    slice_head(n = top_n) %>%
    pull(Category)
  
  # 处理数据，将非top_n的分类归为"Other"
  data <- data %>%
    group_by(Sample) %>%
    mutate(Category = if_else(Category %in% category_sums, Category, "Other")) %>%
    group_by(Sample, Category) %>%
    summarise(Value = sum(Value, na.rm = TRUE), .groups = "drop") %>%
    # 确保所有样本的总和为100%
    group_by(Sample) %>%
    mutate(Percentage = Value / sum(Value, na.rm = TRUE) * 100) %>%
    ungroup()
  
  # 设置颜色
  n_categories <- length(unique(data$Category))
  col5 <- colorRampPalette((pal_npg("nrc")(9)))(n_categories)
  
  ################################ 冲积图 ##################################
  
  # alluvial_plot <- ggplot(data = data, 
  #                         aes(x = Sample, y = Percentage, 
  #                             alluvium = Category, 
  #                             stratum = Category)) +
  #   geom_alluvium(aes(fill = Category), 
  #                 color = NA, alpha = 0.5, 
  #                 decreasing = FALSE, width = 1/2) +
  #   geom_stratum(aes(fill = Category), 
  #                color = NA, decreasing = FALSE, width = 1/2) +
  #   scale_fill_manual(values = rev(col5)) +
  #   labs(y = "Relative abundance (%)", x = "", fill = type) +
  #   theme_minimal() +
  #   theme(
  #     panel.grid.major = element_blank(),
  #     panel.grid.minor = element_blank(),
  #     axis.text.x = element_text(angle = 45, hjust = 1),
  #     axis.ticks.y = element_line(color = "black"),
  #     axis.line.y = element_line(color = "black"),
  #     strip.background = element_rect(fill = "grey95", colour = "grey95")
  #   )
  # 
  # ggsave(file.path(output_folder, paste0(type, "_alluvial_plot.pdf")), 
  #        plot = alluvial_plot, width = 11, height = 8, units = "in")
  # 
  ################################ 堆叠柱状图 ##################################
  
  # 添加分组变量，每50个样本一组
  data <- data %>%
    mutate(Plot_Group = ceiling(as.numeric(factor(Sample)) / 50))  # 每50个样本分为一组
  
  # 计算分组数量
  n_groups <- length(unique(data$Plot_Group))
  # 根据分组数量计算合适的高度：每组大约3英寸，最小6英寸
  plot_height <- max(6, n_groups * 3)
  
  bar_plot <- ggplot(data, aes(x = Sample, y = Percentage, fill = Category)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = rev(col5)) +
    facet_wrap(~ Plot_Group, scales = "free_x", ncol = 1) +  # 每50个样本一组，垂直排列
    labs(y = "Relative abundance (%)", x = "", fill = type) +
    theme_minimal() +
    theme(
      text = element_text(family = main_font),  # 设置全局字体
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.ticks.y = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      strip.background = element_rect(fill = "grey95", colour = "grey95"),
      strip.text = element_blank(),  # 隐藏分组标签
      panel.spacing = unit(1, "lines"),
      plot.title = element_text(size = 14),  # 主标题
      axis.title = element_text(size = 12),   # 坐标轴标题
      axis.text = element_text(size = 10),    # 坐标轴刻度文本
      legend.title = element_text(size = 10), # 图例标题
      legend.text = element_text(size = 10)   # 图例文本
    )
  
  # 使用动态计算的高度保存图片
  ggsave(file.path(output_folder, paste0(type, "_barplot.pdf")), 
         plot = bar_plot, 
         width = 12,
         height = plot_height,  # 使用动态计算的高度
         units = "in")
  
  ################################ 分组图 ##################################
  if (file.exists(group_file)) {
    # 读取分组信息
    data_group <- tryCatch({
      group_data <- read.csv(group_file, sep = '\t', header = TRUE, 
                             check.names = FALSE, 
                             col.names = c("Sample", "Group"),
                             stringsAsFactors = FALSE)
      # 清理分组数据
      group_data <- group_data[!is.na(group_data$Sample) & !is.na(group_data$Group), ]
      group_data$Sample <- as.character(group_data$Sample)
      group_data$Group <- as.character(group_data$Group)
      group_data
    }, error = function(e) {
      stop("读取分组文件时出错: ", e$message)
    })
    
    # 合并数据
    data_with_groups <- data %>%
      inner_join(data_group, by = "Sample")
    
    # 计算每个组的样本数量
    group_sizes <- data_with_groups %>%
      group_by(Group) %>%
      summarise(n_samples = n_distinct(Sample)) %>%
      arrange(desc(n_samples))  # 按样本数量降序排列
    
    # 计算合适的每行分组数
    max_samples_per_group <- max(group_sizes$n_samples)
    # 假设每个样本在图中占用约0.15英寸宽度，考虑到图例和边距
    groups_per_row <- floor(11 / (max_samples_per_group * 0.15))  # 11是图片总宽度
    groups_per_row <- max(1, min(groups_per_row, 4))  # 限制每行最少1个，最多4个分组
    
    # 计算所需的总行数
    n_groups <- length(unique(data_with_groups$Group))
    n_rows <- ceiling(n_groups / groups_per_row)
    
    # 计算合适的图片高度（每行大约需要3英寸）
    plot_height <- max(8, n_rows * 3)
    
    # 计算x轴标签的字体大小
    # 基准：如果每组少于20个样本，使用8pt字体
    # 样本数量越多，字体越小，但不小于4pt
    base_font_size <- 8
    max_samples <- max(group_sizes$n_samples)
    font_size <- max(8, min(12, base_font_size * (20 / max_samples)))
    
    # 绘制分组堆叠柱状图
    grouped_bar_plot <- ggplot(data_with_groups, 
                               aes(x = Sample, y = Percentage, fill = Category)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = rev(col5)) +
      labs(y = "Relative abundance (%)", x = "", fill = type) +
      theme_minimal() +
      theme(
        text = element_text(family = main_font),  # 设置全局字体
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        strip.background = element_rect(fill = "grey95", colour = "grey95"),
        strip.text = element_text(size = 12, face = "bold"),
        panel.spacing = unit(1, "lines"),
        plot.title = element_text(size = 14),  # 主标题
        axis.title = element_text(size = 12),   # 坐标轴标题
        axis.text = element_text(size = 10),    # 坐标轴刻度文本
        legend.title = element_text(size = 10), # 图例标题
        legend.text = element_text(size = 10),   # 图例文本
        axis.text.x = element_text(
          angle = 45, 
          hjust = 1,
          size = font_size  # 使用动态计算的字体大小
        ),
      ) +
      facet_wrap(~ Group, scales = "free_x", 
                 ncol = groups_per_row)  # 使用计算出的每行分组数
    
    # 保存图片，使用计算出的高度
    ggsave(file.path(output_folder, paste0(type, "_grouped_barplot.pdf")), 
           plot = grouped_bar_plot, 
           width = 11, 
           height = plot_height, 
           units = "in")
    
    
    ################################ 分组平均堆叠柱状图 ##################################
    
    # 计算每个组的平均值
    group_mean_data <- data_with_groups %>%
      group_by(Group, Category) %>%
      summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
      group_by(Group) %>%
      mutate(Percentage = Value / sum(Value, na.rm = TRUE) * 100)
    
    # 绘制分组平均堆叠柱状图
    group_mean_plot <- ggplot(group_mean_data, 
                              aes(x = Group, y = Percentage, fill = Category)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = rev(col5)) +
      labs(y = "Relative abundance (%)", x = "", fill = type) +
      theme_minimal() +
      theme(
        text = element_text(family = main_font),  # 设置全局字体
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.y = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        strip.background = element_rect(fill = "grey95", colour = "grey95"),
        plot.title = element_text(size = 14),  # 主标题
        axis.title = element_text(size = 12),   # 坐标轴标题
        axis.text = element_text(size = 10),    # 坐标轴刻度文本
        legend.title = element_text(size = 10), # 图例标题
        legend.text = element_text(size = 10)   # 图例文本
      )
    
    ggsave(file.path(output_folder, paste0(type, "_grouped_mean_barplot.pdf")), 
           plot = group_mean_plot, width = 11, height = 8, units = "in")
  }
}

args <- commandArgs(trailingOnly = TRUE)
type <- args[1]
top_n <- as.numeric(args[2])

input_file <- args[3]
group_file <- args[4]
output_folder <- args[5]


# type <- "Kingdom"
# top_n <- 20
# 
# input_file <- "E:/XuYunProject/R绘图项目/冲积图/单样本_input/Kingdom_abundance1.xls"
# group_file <- "E:/XuYunProject/R绘图项目/冲积图/单样本_input/single_Sample_Group_time.csv"
# output_folder <- "E:/XuYunProject/R绘图项目/冲积图/单样本_output_新罗马/"

# type <- "drug"
# top_n <- 20
# 
# input_file <- "E:/XuYunProject/R绘图项目/冲积图/test_small_input/card_drug_class_abundance.xls"
# group_file <- "E:/XuYunProject/R绘图项目/冲积图/test_small_input/sampleGroup.tsv"
# output_folder <- "E:/XuYunProject/R绘图项目/冲积图/test_small_output/"

# 定义需要绘制的分类级别
plot_abundance(type, input_file, group_file, top_n, output_folder)