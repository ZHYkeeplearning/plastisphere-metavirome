library(ggtree)
library(ape)
library(dplyr)
library(ggplot2)
library(stringr)
library(extrafont)
loadfonts(device = "pdf")

# 函数：将数据转换为Newick格式
convert_to_newick <- function(data) {
  data$parent <- as.character(data$parent)
  data$type <- as.character(data$type)
  
  build_subtree <- function(node, data) {
    children <- data$type[data$parent == node]
    if (length(children) == 0) {
      return(node)
    } else {
      subtrees <- sapply(children, build_subtree, data = data)
      return(paste0("(", paste(subtrees, collapse = ","), ")", node))
    }
  }
  
  root <- "root"
  newick <- build_subtree(root, data)
  newick <- paste0(newick, ";")
  return(newick)
}

# 函数：递归计算每个节点的层级
calculate_depth <- function(node, tree) {
  if (node <= length(tree$tip.label)) {
    return(0)
  }
  
  children <- tree$edge[tree$edge[, 1] == node, 2]
  if (length(children) == 0) {
    return(0)
  } else {
    depths <- sapply(children, calculate_depth, tree = tree)
    return(max(depths) + 1)
  }
}

# 读取命令行参数
args <- commandArgs(trailingOnly = TRUE)
intpu_type <- args[1]
input_file <- args[2]
input_abundance_file <- args[3]
output_dir <- args[4]

# intpu_type <- ifelse(length(args) >= 1, args[1], "Class")
# output_dir <- ifelse(length(args) >= 4, args[4], "E:/XuYunProject/R绘图项目/树状图/output_6_sample/")
# input_file <- ifelse(length(args) >= 2, args[2], "E:/XuYunProject/R绘图项目/树状图/input_6_sample/Class_abundance.xls")
# input_abundance_file <- ifelse(length(args) >= 3, args[3], "E:/XuYunProject/R绘图项目/树状图/input_6_sample/taxonomy_abundance.xls")

# 读取数据
data_abundance <- read.csv(input_file, sep="\t", header=TRUE, check.names = FALSE, row.names = 1, stringsAsFactors = FALSE)
taxonomy_data <- read.csv(input_abundance_file, sep="\t", header=TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# 获取样本的名称
sample_names <- colnames(data_abundance)

# 创建一个空列表来存储每个样本的类型数组
list_of_types <- list()

# 遍历所有样本，为每个样本创建一个类型数组
for (sample in sample_names) {
  current_sample_data <- data.frame(Value = data_abundance[, sample])
  rownames(current_sample_data) <- rownames(data_abundance)
  
  filtered_data <- current_sample_data %>%
    mutate(Value = round(Value, digits = 2)) %>%
    filter(Value >= 0.01) %>%
    arrange(desc(Value)) %>%
    head(n = 30)
  
  types <- rownames(filtered_data)
  list_of_types[[sample]] <- types
}

# 处理每个样本
for (sample in names(list_of_types)) {
  print(paste("Processing sample:", sample))
  input_types <- list_of_types[[sample]]
  
  relationship_df <- data.frame(child = character(), parent = character(), stringsAsFactors = FALSE)
  added_relationships <- character()
  
  for (type in input_types) {
    for (i in 1:nrow(taxonomy_data)) {
      current_taxonomy <- taxonomy_data[i, 1]
      
      # 使用固定字符串匹配而不是正则表达式
      pattern <- paste0("|", type, "|")
      end_pattern <- paste0("|", type)
      
      if (grepl(pattern, current_taxonomy, fixed = TRUE) || 
          endsWith(current_taxonomy, end_pattern)) {
        # 使用fixed=TRUE进行分割
        split_taxonomy <- strsplit(current_taxonomy, "|", fixed = TRUE)[[1]]
        split_taxonomy <- split_taxonomy[split_taxonomy != ""]
        
        for (j in 1:(length(split_taxonomy) - 1)) {
          child <- split_taxonomy[j + 1]
          parent <- split_taxonomy[j]
          relationship_id <- paste(parent, child, sep = "->")
          if (!(relationship_id %in% added_relationships)) {
            relationship_df <- rbind(relationship_df, data.frame(child = child, parent = parent))
            added_relationships <- c(added_relationships, relationship_id)
          }
        }
        break
      }
    }
  }
  
  relationship_df$percentage <- 0
  
  # 计算总和
  sample_total <- taxonomy_data %>%
    summarise(total = sum(!!sym(sample), na.rm = TRUE)) %>%
    pull(total)
  
  # 为每个子节点计算百分比
  relationship_df <- relationship_df %>%
    rowwise() %>%
    mutate(percentage = {
      child_pattern <- paste0("|", child, "|")
      child_end_pattern <- paste0("|", child)
      
      matches <- taxonomy_data %>%
        filter(grepl(child_pattern, taxonomy_data[[1]], fixed = TRUE) |
                 endsWith(taxonomy_data[[1]], child_end_pattern))
      
      sum(matches[[sample]], na.rm = TRUE) / sample_total * 100
    })
  
  # 解除行绑定
  relationship_df <- ungroup(relationship_df)
  
  # 处理父节点
  for (parent_type in unique(relationship_df$parent)) {
    if (grepl("^k__", parent_type)) {
      if (!(parent_type %in% relationship_df$child)) {
        parent_matches <- taxonomy_data %>%
          filter(str_detect(taxonomy_data[[1]], fixed(parent_type)))
        
        total_percentage <- sum(parent_matches[[sample]], na.rm = TRUE) / sample_total * 100
        
        new_row <- data.frame(child = parent_type, parent = "root", percentage = total_percentage)
        relationship_df <- rbind(relationship_df, new_row)
      }
    }
  }
  
  names(relationship_df)[names(relationship_df) == "child"] <- "type"
  relationship_df$percentage <- round(relationship_df$percentage, digits = 2)
  
  colors <- c("#FF6633", "#FFB399", "#FF33FF", "#FFFF99", "#00B3E6", "#E6B333")
  relationship_df$color <- sample(colors, nrow(relationship_df), replace = TRUE)
  
  newick_str <- convert_to_newick(relationship_df)
  tree <- read.tree(text = newick_str)
  
  node_depths <- sapply(1:(length(tree$tip.label) + tree$Nnode), calculate_depth, tree = tree)
  node_names <- c(tree$tip.label, tree$node.label)
  depth_df <- data.frame(node = node_names, depth = node_depths)
  
  node_colors <- colors[(node_depths %% length(colors)) + 1]
  
  labels <- relationship_df %>%
    mutate(label = paste0(type, "\n", percentage, "%"))
  
  node_mapping <- setNames(relationship_df$percentage, relationship_df$type)
  
  node_color_mapping <- data.frame(node = node_names, color = node_colors)
  
  p <- ggtree(tree) +
    geom_tiplab(aes(label = relationship_df$type[match(label, relationship_df$type)]), color = "black", vjust = -0.5, hjust = 1, nudge_x = -0.1, size = 2.5) +
    geom_text2(aes(subset = !isTip & label != "root", label = relationship_df$type[match(label, relationship_df$type)]), color = "black", vjust = -0.5, hjust = 1, nudge_x = -0.1, size = 2.5) +
    geom_text2(aes(subset = !is.na(node_mapping[label]) & label != "root", label = paste0(node_mapping[label], "%")), color = "black", vjust = 1.5, hjust = 1, nudge_x = -0.1, size = 2.5) +
    geom_point2(aes(subset = !is.na(node_mapping[label]) & label != "root", size = node_mapping[label], fill = node_color_mapping$color[match(label, node_color_mapping$node)]), shape = 21) +
    scale_size_continuous(range = c(1, 10), guide = "none") +
    scale_fill_manual(values = unique(node_color_mapping$color), guide = "none")
  
  ggsave(paste0(output_dir, sample, '.', intpu_type,'_tree.pdf'), plot = p, width = 10, height = 20, units = "in", family = "Times New Roman")
}

