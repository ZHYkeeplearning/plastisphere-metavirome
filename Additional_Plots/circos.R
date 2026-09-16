# 加载包
library(circlize)
library(tidyr)
library(tibble)
library(dplyr)
library(extrafont)
loadfonts(device = "pdf")

# 处理命令行参数
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
group_file <- args[2]
output_pdf <- args[3]
output_xls <- args[4]

# # （测试时取消注释这些行）
# input_file <- "E:/XuYunProject/R绘图项目/Circos/input/Kindom_abundance.xls"
# group_file <- "E:/XuYunProject/R绘图项目/Circos/input/group.xls"
# output_pdf <- "E:/XuYunProject/R绘图项目/Circos/output/Kindom_circos.pdf"
# output_xls <- "E:/XuYunProject/R绘图项目/Circos/output/Kindom_circos.xls"

# 读取数据
data_abundance <- read.csv(input_file, sep = "\t", header = TRUE, check.names = FALSE, row.names = 1)

# 筛选前10个高丰度物种
data_abundance$sum <- rowSums(data_abundance)
top_data <- head(data_abundance[order(data_abundance$sum, decreasing = TRUE), ], 10)
top_data <- top_data[, -ncol(top_data)]
top_data <- rownames_to_column(top_data, var = "Taxonomy")

# 读取分组信息
data_group <- read.csv(group_file, sep = '\t', header = TRUE, 
                       check.names = FALSE, 
                       col.names = c("Sample", "Group")) %>%
  mutate(across(c(Sample, Group), as.character))

# 数据长格式转换
data_combined <- top_data %>%
  pivot_longer(cols = -Taxonomy, names_to = "Sample", values_to = "Abundance") %>%
  left_join(data_group, by = "Sample") %>%
  rename(group = Group)

# 输出中间结果
write.table(data_combined, output_xls, sep = '\t', quote = FALSE, row.names = FALSE)

# 准备绘图矩阵
wide_data <- data_combined %>%
  select(-Sample) %>%
  group_by(Taxonomy, group) %>%
  summarise(Abundance = sum(Abundance), .groups = 'drop') %>%
  pivot_wider(names_from = group, values_from = Abundance, values_fill = list(Abundance = 0))

mat_wide_data <- as.matrix(wide_data[, -1])
rownames(mat_wide_data) <- wide_data$Taxonomy

# 设置PDF输出（关键字体设置）
pdf(output_pdf, width = 12, height = 12, family = "Times New Roman")  # 设置PDF字体
par(mar = c(1, 1, 1, 1), oma = c(1, 1, 1, 1), family = "Times New Roman")  # 设置全局字体

# 调整绘图参数
circos.par(canvas.xlim = c(-1.5, 1.5), 
           canvas.ylim = c(-1.5, 1.5),
           gap.after = structure(rep(5, nrow(mat_wide_data) + ncol(mat_wide_data)),
                                 names = c(rownames(mat_wide_data), colnames(mat_wide_data))))

# 绘制弦图
chordDiagram(mat_wide_data,
             annotationTrack = "grid",
             transparency = 0.5,
             preAllocateTracks = list(track.height = 0.2),
             reduce = 0)

# 添加坐标轴
for(si in get.all.sector.index()) {
  circos.axis(h = "top", 
              labels.cex = 0.6,
              sector.index = si,
              track.index = 2,
              labels.font = 1)  # 设置坐标轴字体为常规
}

# 添加扇区标签
circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), 
                mean(ylim), 
                sector.name,
                niceFacing = TRUE,
                facing = "reverse.clockwise",
                adj = c(1, 0.5),
                cex = 1,
                font = 1)  # 设置标签字体为常规
  },
  bg.border = NA
)

# 清理绘图设备
circos.clear()
dev.off()

print(paste("结果已保存至:", output_pdf))