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


draw_pca = function(inputFile, groupFile, outputFile, addEllipse) {
  tryCatch({
    # 读取数据部分保持不变
    pca_data <- fread(inputFile, header=TRUE, sep="\t", 
                      check.names=FALSE, 
                      quote="", 
                      data.table=FALSE)
    
    rownames(pca_data) <- pca_data[[1]]
    pca_data <- pca_data[, -1]
    
    pca_data <- t(pca_data)
    pca_data <- data.frame(pca_data, check.names=FALSE)
    
    group_list <- fread(groupFile, header = TRUE, sep = "\t",
                        check.names = FALSE, 
                        col.names = c("Sample", "Group"))
    
    pca_data$Group <- as.factor(group_list$Group)
    
    # 执行PCA
    pca <- PCA(pca_data[, -ncol(pca_data), drop = FALSE], 
               graph = FALSE)
    
    # 创建绘图数据框
    pca_df <- data.frame(
      Ind = rownames(pca$ind$coord),
      Dim1 = pca$ind$coord[, 1],
      Dim2 = pca$ind$coord[, 2],
      Group = pca_data$Group
    )
    
    # 扩展形状和颜色集合
    custom_shapes <- c(
      16, # 实心圆
      17, # 实心三角形
      15, # 实心方块
      18, # 实心菱形
      8,  # 星形
      3,  # 加号
      4,  # 叉号
      5,  # 菱形
      6,  # 倒三角
      7,  # 方块带x
      9,  # 菱形带+
      10, # 三角形带点
      11, # 倒三角带点
      12, # 方块带点
      13  # 圆形带点
    )
    
    
    custom_colors <- c(
      "#E41A1C", # 红色
      "#377EB8", # 蓝色
      "#4DAF4A", # 绿色
      "#984EA3", # 紫色
      "#FF7F00", # 橙色
      "#FF5000", # 橘色
      "#A65628", # 棕色
      "#F781BF", # 粉色
      "#999999", # 灰色
      "#66C2A5", # 青色
      "#FC8D62", # 橙红
      "#8DA0CB", # 淡蓝
      "#E78AC3", # 深粉
      "#A6D854", # 黄绿
      "#FFD92F"  # 金黄
    ) 
    n_groups <- length(unique(pca_df$Group))
    
    # 确保有足够的颜色和形状
    if(n_groups > length(custom_colors)) {
      custom_colors <- rep(custom_colors, ceiling(n_groups/length(custom_colors)))[1:n_groups]
    }
    if(n_groups > length(custom_shapes)) {
      custom_shapes <- rep(custom_shapes, ceiling(n_groups/length(custom_shapes)))[1:n_groups]
    }
    
    # 创建基础图形
    p <- ggplot(pca_df, aes(x = Dim1, y = Dim2, color = Group, shape = Group)) +
      geom_point(size = 4) +  # 增大点的大小
      scale_shape_manual(values = custom_shapes[1:n_groups]) +
      scale_color_manual(values = custom_colors[1:n_groups]) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "#333333") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#333333") +
      xlab(paste("PC1 (", round(pca$eig[1, 2], 1), "%)", sep = "")) +
      ylab(paste("PC2 (", round(pca$eig[2, 2], 1), "%)", sep = "")) +
      theme_bw() +
      theme(
        text = element_text(family = main_font),  # 设置全局字体
        panel.grid = element_blank(),
        legend.position = "right",
        # 调整图例布局
        legend.key.size = unit(1.2, "lines"),
        legend.spacing.y = unit(0.5, "lines"),
        legend.box.spacing = unit(1, "lines"),
        plot.title = element_text(size = 14),  # 主标题
        axis.title = element_text(size = 12),   # 坐标轴标题
        axis.text = element_text(size = 10),    # 坐标轴刻度文本
        legend.title = element_text(size = 10), # 图例标题
        legend.text = element_text(size = 10)   # 图例文本
      )
    
    # 如果组别较多，调整图例为多列显示
    if(n_groups > 8) {
      p <- p + guides(
        color = guide_legend(ncol = 2),
        shape = guide_legend(ncol = 2)
      )
    }
    
    if (addEllipse) {
      p <- p + stat_ellipse(type = "t", level = 0.95, 
                            aes(color = Group), alpha = 0.2)
    }
    
    # 保存图片，如果组别多，适当增加宽度
    width <- if(n_groups > 8) 12 else 10
    ggsave(outputFile, plot = p, width = width, height = 8, dpi = 300)
    
    return(p)
    
  }, error = function(e) {
    print(paste("Error occurred:", e$message))
    return(NULL)
  })
}


args <- commandArgs(trailingOnly = TRUE)

inputFile <- args[1]
groupFile <- args[2]
outputFile <- args[3]
addEllipse <- as.numeric(args[4])
# 
# 使用示例
# inputFile <- "E:/XuYunProject/R绘图项目/PCA/input_error/Species_abundance.xls"
# groupFile <- "E:/XuYunProject/R绘图项目/PCA/input_error/FLOWHUB_202400070_1_SampleGroup.csv"
# outputFile <- "E:/XuYunProject/R绘图项目/PCA/output_error/dist_PCA1.pdf"
# addEllipse <- 0

# inputFile <- "E:/XuYunProject/R绘图项目/PCA/test_small_input/card_drug_class_abundance.xls"
# groupFile <- "E:/XuYunProject/R绘图项目/PCA/test_small_input/sampleGroup.tsv"
# outputFile <- "E:/XuYunProject/R绘图项目/PCA/test_small_output/dist_PCA1.pdf"
# addEllipse <- 0

result <- draw_pca(inputFile, groupFile, outputFile, addEllipse)

if (is.null(result)) {
  print("PCA execution or plotting failed. Please check input data and error messages.")
} else {
  print("PCA execution and plotting completed successfully.")
}
