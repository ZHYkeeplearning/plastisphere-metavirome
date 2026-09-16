library(pheatmap)
library(ggplot2)
library(extrafont)
loadfonts(device = "pdf")
# calcMethod 可选默认为pearson皮尔逊相关系数,kendall肯德尔秩相关系数, spearman 斯皮尔曼系数

#
args <- commandArgs(trailingOnly = TRUE)
dataFile = args[1]
groupFile = args[2]
outputFile = args[3]

# # # 设置参数
# dataFile = "E:/XuYunProject/R绘图项目/热图/input_单样本/Class_count.xls"
# groupFile = "E:/XuYunProject/R绘图项目/热图/input_单样本/SampleGroup.txt"
# outputFile = "E:/XuYunProject/R绘图项目/热图/output_单样本/Class_count_grouped_hetamap.pdf"

# # 设置参数
# dataFile = "E:/XuYunProject/R绘图项目/热图/多样本_input/Kingdom_abundance.xls"
# groupFile = "E:/XuYunProject/R绘图项目/热图/多样本_input/180_Sample_Group_time.csv"
# outputFile = "E:/XuYunProject/R绘图项目/热图/多样本_output/Kingdom_count_grouped_hetamap.pdf"

class_group = ""
type = "sample_row"
showValue = FALSE
cluster_row = FALSE
cluster_col = FALSE
calcMethod = "pearson"

# 读取输入数据
df = read.csv(dataFile, sep = "\t", header = TRUE, check.names = FALSE, row.names = 1)

# 删除 Detail 列（如果存在）
if ("Detail" %in% colnames(df)) {
  df <- df[, !colnames(df) %in% "Detail"]
}

# 保存行名
row_names <- rownames(df)

# 转换为数值类型，但保持行名
df <- as.data.frame(apply(df, 2, function(x) as.numeric(as.character(x))))
rownames(df) <- row_names

# 如果需要计算相关性
if (type == "sample") {
  if (calcMethod == "pearson") {
    pheatmap_data <- cor(df)
  } else if (calcMethod == "kendall") {
    pheatmap_data <- cor(df, method = "kendall")
  } else if (calcMethod == "spearman") {
    pheatmap_data <- cor(df, method = "spearman")
  }
} else {
  pheatmap_data <- df
}

# 计算每一行的总和
pheatmap_data$sum <- rowSums(pheatmap_data, na.rm = TRUE)

# 对总和进行排序，并筛选出总和不为0的类别
top_data <- pheatmap_data[pheatmap_data$sum > 0, ]
top_data <- top_data[order(top_data$sum, decreasing = TRUE), ]

# 删除总和列
top_data <- top_data[, -ncol(top_data), drop = FALSE]

# 如果类别数量超过10个，只保留前10个
if (nrow(top_data) > 10) {
  top_data <- top_data[1:10, ]
}

# 转置数据以确保样本名作为行名
plot_data <- t(top_data)

# 读取分组信息（如果提供）
if (groupFile != "") {
  groupList = read.csv(groupFile, sep = "\t", header = TRUE, row.names = 1)
} else {
  groupList = NULL
}

# 检查分组信息是否与样本对应
if (!is.null(groupList)) {
  groupList <- groupList[rownames(plot_data), , drop = FALSE]
}

# 确保分组顺序与图例一致
if (!is.null(groupList)) {
  for (colname in colnames(groupList)) {
    groupList[[colname]] <- factor(groupList[[colname]], levels = unique(groupList[[colname]]))
  }
}

# 读取类别信息（如果提供）
if (class_group != "") {
  classGroup = read.csv(class_group, sep = "\t", header = TRUE, row.names = 1)
  col_class <- classGroup[colnames(plot_data), , drop = FALSE]
} else {
  col_class = NULL
}

# 计算动态宽度和高度
base_height = 10
sample_count = nrow(plot_data)
height_factor = ceiling(sample_count / 10) * 2
plot_height = base_height + height_factor
plot_width = 12

# 创建自定义注释颜色
ann_colors = list()
if (!is.null(groupList)) {
  for (colname in colnames(groupList)) {
    if (is.factor(groupList[[colname]]) || is.character(groupList[[colname]])) {
      unique_values <- unique(groupList[[colname]])
      colors <- rainbow(length(unique_values))
      ann_colors[[colname]] <- setNames(colors, unique_values)
    }
  }
}
# 绘制热图时添加字体参数 (关键修改)
p <- pheatmap::pheatmap(plot_data, 
                        display_numbers = showValue, 
                        annotation_row = groupList,
                        annotation_col = col_class,  
                        cluster_cols = cluster_row,  
                        cluster_rows = cluster_col,
                        fontfamily = "Times New Roman",  # 设置字体家族
                        fontsize_row = ifelse(sample_count > 50, 6,
                                              ifelse(sample_count > 30, 8,
                                                     ifelse(sample_count > 20, 10, 12))),
                        fontsize_col = 10,  # 显式设置列名字体大小
                        annotation_legend = TRUE,
                        annotation_names_row = TRUE,
                        annotation_names_col = TRUE,
                        show_rownames = TRUE,
                        show_colnames = TRUE,
                        angle_col = 270,
                        annotation_colors = ann_colors,
                        gaps_col = 0,
                        # 以下为新增的图例字体设置
                        legend_labels = element_text(family = "Times New Roman"),
                        legend_title = element_text(family = "Times New Roman"))

# 保存PDF时使用字体嵌入 (关键修改)
pdf(outputFile, width = plot_width, height = plot_height,
    family = "Times New Roman") # 确保使用正确的字体名称

# 打印图形对象
print(p)

# 关闭图形设备
dev.off()