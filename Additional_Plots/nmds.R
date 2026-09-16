library(data.table)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(ape)
library(vegan)
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


draw_nmds = function(inputFile, groupFile, outputFile, distanceType, addEllipse) {
  # 自定义形状和颜色
  custom_shapes <- c(
    16, 17, 15, 18, 8, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13
  )
  
  custom_colors <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
    "#FF5000", "#A65628", "#F781BF", "#999999", "#66C2A5",
    "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F"
  )
  
  # 读取数据并进行错误处理
  pca_data <- tryCatch({
    read.table(inputFile, header=TRUE, sep="\t", row.names=1, check.names=FALSE)
  }, error = function(e) {
    print(paste("Error reading input file:", e$message))
    return(NULL)
  })
  
  if (is.null(pca_data)) return(NULL)
  
  # NA值处理
  print("Checking NA values before transpose...")
  na_count <- sum(is.na(pca_data))
  print(paste("Number of NA values:", na_count))
  
  if (na_count > 0) {
    print("Original dimensions:")
    print(dim(pca_data))
    pca_data <- na.omit(pca_data)
    print("Dimensions after NA removal:")
    print(dim(pca_data))
  }
  
  # 数据转换和处理
  pca_data <- t(pca_data)
  pca_data <- data.frame(pca_data, check.names=FALSE)
  pca_data <- apply(pca_data, 2, function(x) as.numeric(as.character(x)))
  pca_data <- data.frame(pca_data, check.names=FALSE)
  
  # 读取分组信息
  group_list <- tryCatch({
    read.table(groupFile, header = TRUE, sep = "\t", 
               check.names = FALSE, col.names = c("Sample", "Group"))
  }, error = function(e) {
    print(paste("Error reading group file:", e$message))
    return(NULL)
  })
  
  if (is.null(group_list)) return(NULL)
  
  # 数据维度检查
  if (nrow(pca_data) != nrow(group_list)) {
    print("Row number mismatch between data and group list")
    print(paste("Data rows:", nrow(pca_data)))
    print(paste("Group list rows:", nrow(group_list)))
    return(NULL)
  }
  
  pca_data$Group <- as.factor(group_list$Group)
  n_groups <- length(unique(pca_data$Group))
  
  # 样本数量检查
  if (nrow(pca_data) <= 3 || ncol(pca_data) <= 3) {
    print("Insufficient samples or variables for NMDS")
    return(NULL)
  }
  
  # NMDS计算和绘图
  tryCatch({
    # 距离矩阵计算
    dist_matrix <- switch(distanceType,
                          "bray" = vegdist(pca_data[, -ncol(pca_data)], method = "bray"),
                          "euclidean" = dist(pca_data[, -ncol(pca_data)], method = "euclidean"),
                          "jaccard" = vegdist(pca_data[, -ncol(pca_data)], method = "jaccard"),
                          stop("Unsupported distance type"))
    
    # NMDS计算
    nmds_result <- metaMDS(dist_matrix, k = 2, trymax = 100)
    stress_value = round(nmds_result$stress, 4)
    
    # 准备绘图数据
    nmds_df <- as.data.frame(nmds_result$points[, 1:2])
    colnames(nmds_df) <- c("NMDS1", "NMDS2")
    nmds_df$Group <- pca_data$Group
    
    # 确保有足够的颜色和形状
    if(n_groups > length(custom_colors)) {
      custom_colors <- rep(custom_colors, ceiling(n_groups/length(custom_colors)))[1:n_groups]
    }
    if(n_groups > length(custom_shapes)) {
      custom_shapes <- rep(custom_shapes, ceiling(n_groups/length(custom_shapes)))[1:n_groups]
    }
    
    # 创建基础图形
    p <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Group, shape = Group)) +
      geom_point(size = 4) +
      scale_shape_manual(values = custom_shapes[1:n_groups]) +
      scale_color_manual(values = custom_colors[1:n_groups]) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "#333333") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#333333") +
      annotate("text", x = max(nmds_df$NMDS1), y = max(nmds_df$NMDS2),
               label = paste("Stress:", stress_value),
               hjust = 1.1, vjust = 1.1, size = 4) +
      theme_bw() +
      theme(
        text = element_text(family = main_font),  # 设置全局字体
        panel.grid = element_blank(),
        legend.position = "right",
        legend.key.size = unit(1.2, "lines"),
        plot.title = element_text(size = 14),  # 主标题
        axis.title = element_text(size = 12),   # 坐标轴标题
        axis.text = element_text(size = 10),    # 坐标轴刻度文本
        legend.title = element_text(size = 10), # 图例标题
        legend.text = element_text(size = 10)   # 图例文本
      )
    
    # 多组别时调整图例
    if(n_groups > 8) {
      p <- p + guides(
        color = guide_legend(ncol = 2),
        shape = guide_legend(ncol = 2)
      )
    }
    
    # 添加置信椭圆
    if (addEllipse) {
      p <- p + stat_ellipse(type = "t", level = 0.95, 
                            aes(color = Group), alpha = 0.2)
    }
    
    # 保存图片
    width <- if(n_groups > 8) 12 else 10
    ggsave(outputFile, plot = p, width = width, height = 8, dpi = 300)
    
    return(p)
    
  }, error = function(e) {
    print(paste("NMDS calculation error:", e$message))
    return(NULL)
  })
}

# 命令行参数处理
args <- commandArgs(trailingOnly = TRUE)
distanceType <- args[1]
inputFile <- args[2]
groupFile <- args[3]
outputFile <- args[4]
addEllipse <- as.numeric(args[5])

# 执行NMDS分析并处理结果
result <- draw_nmds(inputFile, groupFile, outputFile, distanceType, addEllipse)

if (is.null(result)) {
  print("NMDS analysis failed. Please check error messages above.")
} else {
  print("NMDS analysis completed successfully.")
}
