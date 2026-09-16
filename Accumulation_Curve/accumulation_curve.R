# =============================================================================
# 脚本名称: plot_accumulation_curve_optimized.R
# 功能: 高效地从大型丰度表 (TSV, 特征 x 样本) 计算并绘制物种累积曲线。
#       - 包含文件读取和转置功能。
#       - 【核心优化】通过将丰度矩阵转换为二进制矩阵，极大提升了计算速度。
#       - 包含将最终图形保存为 PDF 的功能。
# =============================================================================


# --- 步骤 0: 安装并加载所有必需的 R 包 ---

# 检查并安装包
if (!requireNamespace("vegan", quietly = TRUE)) install.packages("vegan")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")

# 加载包
library(vegan)
library(ggplot2)
library(dplyr)
library(data.table)


#-------------------------------------------------------------------------------
# 步骤 1: 准备数据 (读取和转置 '特征 x 样本' TSV 文件)
#-------------------------------------------------------------------------------

# --- 1.1 定义一个函数，用于读取 TSV 文件并将其转置 ---
read_and_transpose_abundance <- function(file_path) {
  cat(sprintf("正在读取并转置文件: %s\n", basename(file_path)))
  
  # 使用 fread 高效读取 TSV 文件
  abundance_dt <- fread(file_path, header = TRUE, data.table = FALSE)
  
  # 将第一列 (特征ID) 设置为行名
  rownames(abundance_dt) <- abundance_dt[[1]]
  
  # 删除原始的第一列
  abundance_dt <- abundance_dt[, -1]
  
  # 执行转置操作，得到 '样本 x 特征' 格式
  abundance_transposed <- as.data.frame(t(abundance_dt))
  
  cat("文件读取和转置成功！\n\n")
  return(abundance_transposed)
}


# --- 1.2 使用上面的函数加载您的三个真实数据文件 ---
# !!! 请将下面的文件路径替换为您真实文件的路径 !!!

# 加载 vPCs 丰度表
vpcs_abundance <- read_and_transpose_abundance(
  file_path = "input/vPCs_contig_count.tsv" # <-- 【请修改】
)

# 加载 vOTUs 丰度表
votus_abundance <- read_and_transpose_abundance(
  file_path = "input/vOTU_contig_count.tsv" # <-- 【请修改】
)

# 加载 mOTUs 丰度表
motus_abundance <- read_and_transpose_abundance(
  file_path = "input/mOTU_contig_count.tsv" # <-- 【请修改】
)


#-------------------------------------------------------------------------------
# 【新增优化步骤】: 将丰度矩阵转换为二进制 (存在/不存在) 矩阵
#-------------------------------------------------------------------------------
# 这一步将极大地加速后续的 specaccum 计算

cat("正在将大型丰度矩阵转换为高效的二进制矩阵...\n")

# 对于 vPCs
vpcs_abundance_binary <- (vpcs_abundance > 0) + 0 # 高效的转换方法: 逻辑判断变为 TRUE/FALSE，+0 将其变为 1/0

# 对于 vOTUs
votus_abundance_binary <- (votus_abundance > 0) + 0

# 对于 mOTUs
motus_abundance_binary <- (motus_abundance > 0) + 0

cat("二进制转换完成！\n\n")


#-------------------------------------------------------------------------------
# 步骤 2 (已更新): 使用【二进制矩阵】进行快速计算
#-------------------------------------------------------------------------------
# 我们现在将转换后的二进制矩阵作为输入，这将非常快

print("正在计算vPCs累积曲线 (快速模式)...")
vpcs_accum <- specaccum(vpcs_abundance_binary, method = "random", permutations = 100)

print("正在计算vOTUs累积曲线 (快速模式)...")
votus_accum <- specaccum(votus_abundance_binary, method = "random", permutations = 100)

print("正在计算mOTUs累积曲线 (快速模式)...")
motus_accum <- specaccum(motus_abundance_binary, method = "random", permutations = 100)


#-------------------------------------------------------------------------------
# 步骤 3: 提取计算结果并整合成一个用于 ggplot2 绘图的数据框
#-------------------------------------------------------------------------------
# 定义一个函数来方便地提取数据
prepare_plot_data <- function(accum_object, category_name, n_permutations) {
  data.frame(
    samples = accum_object$sites,
    mean = accum_object$richness,
    sd = accum_object$sd,
    sem = accum_object$sd / sqrt(n_permutations), 
    category = category_name
  )
}

# 设定我们在步骤2中使用的置换次数
PERMUTATIONS_COUNT <- 100

# 提取并合并数据
plot_data <- bind_rows(
  prepare_plot_data(vpcs_accum, "vPCs", PERMUTATIONS_COUNT),
  prepare_plot_data(votus_accum, "vOTUs", PERMUTATIONS_COUNT),
  prepare_plot_data(motus_accum, "mOTUs", PERMUTATIONS_COUNT)
)


#-------------------------------------------------------------------------------
# 步骤 4: 根据图中描述对数据进行缩放

plot_data_scaled <- plot_data %>%
  mutate(
    # 使用 case_when 来根据不同类别应用不同缩放因子
    scaling_factor = case_when(
      category == "mOTUs" ~ 10,  
      category == "vPCs" ~ 10,  
      TRUE ~ 1                 
    ),
    # 对均值和误差应用上面计算出的缩放因子
    mean_scaled = mean / scaling_factor,
    sem_scaled = sem / scaling_factor
  )

#-------------------------------------------------------------------------------
# 步骤 5 (已更新): 使用 ggplot2 绘制带有误差棒的图形
#-------------------------------------------------------------------------------
# 定义颜色和图例标签
custom_colors <- c("vPCs" = "#89D0D3", "vOTUs" = "#F9B98A", "mOTUs" = "#B8B8B8")
custom_labels <- c(
  "vPCs" = "vPCs (×10)",          
  "vOTUs" = "vOTUs",        
  "mOTUs" = "mOTUs (×10)"    
)

final_plot <- ggplot(plot_data_scaled, aes(x = samples, y = mean_scaled, color = category)) +
  
  # 1. 【核心修改】使用 geom_errorbar 替换 geom_ribbon
  #    它会在每个数据点上绘制垂直的误差棒
  geom_errorbar(
    aes(ymin = mean_scaled - sem_scaled, ymax = mean_scaled + sem_scaled),
    width = 0.8,     # 控制误差棒顶部和底部横杠的宽度
    linewidth = 0.5, # 控制误差棒线条的粗细
    alpha = 0.6      # 设置一点透明度，让线条不那么突兀
  ) +
  
  # 2. 绘制均值曲线 (在线条层之下绘制误差棒，所以后画)
  geom_line(linewidth = 1) +
  
  # 3. 绘制数据点 (在最上层，覆盖住线条和误差棒的中心)
  geom_point(size = 1.5, aes(fill = category), shape = 21, stroke = 0.5, color = "white") +
  
  # 4. 自定义颜色和图例
  #    注意：color 用于线条和误差棒，fill 用于数据点的填充
  scale_color_manual(values = custom_colors, labels = custom_labels, name = "") +
  scale_fill_manual(values = custom_colors, labels = custom_labels, name = "") +
  
  # 5. 设置坐标轴标签和标题
  labs(
    x = "# of samples",
    y = "Cumulative number"
  ) +
  
  # 6. 使用一个简洁的主题，并美化细节
  theme_classic() +
  theme(
    legend.position = "top",
    legend.title = element_blank(), # 彻底隐藏图例标题
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black")
  )

# 显示图形
print(final_plot)

#-------------------------------------------------------------------------------
# 步骤 6: 保存图形为 PDF 文件
#-------------------------------------------------------------------------------
# !!! 请修改为您希望保存文件的目录和文件名 !!!
output_directory <- "output"
output_filename <- "accumulation_curve.pdf"

# 检查输出目录是否存在，如果不存在则创建
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

# 构造完整的文件路径
full_path <- file.path(output_directory, output_filename)

# 使用 ggsave 保存文件
ggsave(
  filename = full_path,
  plot = final_plot,
  width = 7,  # 宽度 (英寸)
  height = 6, # 高度 (英寸)
  units = "in",
  device = 'pdf'
)

cat(sprintf("\n图形已成功保存到: %s\n", full_path))


#-------------------------------------------------------------------------------
# 步骤 7: 【新增】保存关键数据表格为 CSV 文件
#-------------------------------------------------------------------------------
cat("\n正在保存分析过程中生成的关键数据表格...\n")

# 1. 定义输出文件名
raw_data_filename <- file.path(output_directory, "accumulation_data_raw.csv")
scaled_data_filename <- file.path(output_directory, "accumulation_data_for_plotting.csv")

# 2. 保存未经缩放的原始计算结果 (plot_data)
#    使用 fwrite 以获得最佳性能
fwrite(plot_data, file = raw_data_filename, row.names = FALSE)
cat(sprintf("  - 原始累积曲线数据已保存到: %s\n", raw_data_filename))


# 3. 保存用于最终绘图的缩放后数据 (plot_data_scaled)
fwrite(plot_data_scaled, file = scaled_data_filename, row.names = FALSE)
cat(sprintf("  - 用于绘图的缩放后数据已保存到: %s\n", scaled_data_filename))


cat("\n所有文件保存完毕！\n")