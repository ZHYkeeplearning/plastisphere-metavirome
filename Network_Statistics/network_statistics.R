# ============================================================================
# 脚本目标: 批量计算所有处理组的网络属性和稳健性，
#           稳健性计算逻辑已更新以匹配示例文件，
#           并按指定格式保存两个最终的汇总文件。
# ============================================================================

# --- 0. 安装和加载必要的 R 包 ---
# install.packages(c("WGCNA", "igraph", "dplyr", "tidyr", "ggClusterNet"))
library(WGCNA)
library(igraph)
library(ggClusterNet)
library(dplyr)
library(tidyr)

# ============================================================================
#          --- 辅助函数: 计算网络脆弱性与【新的】稳健性 ---
# ============================================================================

# (函数 1, 新增!) 单次稳健性模拟 (逻辑来自您的参考脚本)
# 这个函数模拟一次随机移除，并考虑次级灭绝
rand.remov.once <- function(net_matrix, sp_ra, remove_fraction = 0.5) {
  num_to_remove <- round(nrow(net_matrix) * remove_fraction)
  if (num_to_remove == 0) return(1)
  
  # 随机选择要移除的节点
  id_rm <- sample(1:nrow(net_matrix), num_to_remove)
  
  # 创建一个副本进行操作
  net_damaged <- net_matrix
  
  # 切断被移除节点的所有连接
  net_damaged[id_rm, ] <- 0
  net_damaged[, id_rm] <- 0
  
  # 计算每个节点的平均相互作用强度 (丰度加权)
  net_strength <- net_damaged * sp_ra
  sp_mean_interaction <- colMeans(net_strength)
  
  # 找出发生次级灭绝的节点 (包括最初被移除的节点)
  # 注意：这里我们计算的是最终幸存的节点
  id_survive <- which(sp_mean_interaction > 0)
  
  # 计算幸存比例
  remain_percent <- length(id_survive) / nrow(net_matrix)
  
  return(remain_percent)
}

# (函数 2) 计算 Vulnerability (无变化)
calculate_vulnerability <- function(g) {
  if (!is_connected(g)) {
    g <- induced_subgraph(g, V(g)[components(g)$membership == which.max(components(g)$csize)])
  }
  if (gorder(g) <= 1) return(0)
  initial_efficiency <- global_efficiency(g, weights = abs(E(g)$weight))
  if (initial_efficiency == 0) return(0)
  node_names <- V(g)$name
  vulnerability_scores <- numeric(length(node_names))
  for (i in 1:length(node_names)) {
    g_damaged <- delete_vertices(g, node_names[i])
    damaged_efficiency <- global_efficiency(g_damaged, weights = abs(E(g_damaged)$weight))
    vulnerability_scores[i] <- (initial_efficiency - damaged_efficiency) / initial_efficiency
  }
  return(mean(vulnerability_scores)) 
}


# ============================================================================
#                   --- 主分析流程 ---
# ============================================================================

# --- 1. 读取输入文件 ---
dir.create("output", showWarnings = FALSE, recursive = TRUE)
host_table <- read.table(file.path("input", "filtered_host_order_abundance.csv"), header = TRUE, row.names = 1, sep = ",", check.names = FALSE)
virus_table <- read.table(file.path("input", "viral_Family_abundance_tpm.csv"), header = TRUE, row.names = 1, sep = ",", check.names = FALSE)
host_names <- rownames(host_table); virus_names <- rownames(virus_table)
otu_table <- rbind(host_table, virus_table)
metadata <- read.table(file.path("input", "new_sample_group.tsv"), header = TRUE, sep = "\t")

# --- 2. 循环处理所有处理组 ---
cat("--- 开始批量处理所有处理组 ---\n")
all_treatments <- unique(metadata$Treatment)
all_properties_results <- list()
all_robustness_results <- list()

for (current_treatment in all_treatments) {
  cat(paste("\n>>> 正在处理组:", current_treatment, "<<<\n"))
  
  # a. 构建母网络
  target_samples <- metadata$SampleID[metadata$Treatment == current_treatment]
  otu_table_treatment <- otu_table[, target_samples]
  n_replicates <- 2 
  present_counts <- rowSums(otu_table_treatment > 0)
  otu_table_filtered <- otu_table_treatment[present_counts >= n_replicates, ]
  if (nrow(otu_table_filtered) < 2) { cat("警告: 物种过少，跳过此组。\n"); next }
  data_for_cor <- t(otu_table_filtered)
  cor_result <- corAndPvalue(data_for_cor, method = "spearman", use = "pairwise.complete.obs")
  cor_matrix <- cor_result$cor; p_matrix <- cor_result$p
  p_adj_vec <- p.adjust(p_matrix[lower.tri(p_matrix)], method = "fdr")
  p_adj_matrix <- matrix(0, nrow=ncol(p_matrix), ncol=ncol(p_matrix)); p_adj_matrix[lower.tri(p_adj_matrix)] <- p_adj_vec; p_adj_matrix <- p_adj_matrix + t(p_adj_matrix)
  adj_matrix <- ifelse(abs(cor_matrix) > 0.6 & p_adj_matrix < 0.01, cor_matrix, 0)
  diag(adj_matrix) <- 0
  mother_network <- graph_from_adjacency_matrix(adj_matrix, mode="undirected", weighted=TRUE, diag=FALSE)
  mother_network <- delete.vertices(mother_network, which(degree(mother_network) == 0))
  V(mother_network)$type <- ifelse(V(mother_network)$name %in% host_names, "Host", "Virus")
  if (gorder(mother_network) == 0) { cat("警告: 无法构建母网络，跳过此组。\n"); next }
  
  # (*** 关键修改: 为稳健性计算准备数据 ***)
  # 1. 获取母网络的加权邻接矩阵
  mother_adj_matrix <- as_adjacency_matrix(mother_network, attr = "weight", sparse = FALSE)
  # 2. 计算母网络中节点的平均相对丰度
  avg_abundance_all <- rowMeans(otu_table_filtered[V(mother_network)$name, ])
  sp_ra <- avg_abundance_all / sum(avg_abundance_all) # 标准化
  
  # b. 提取子网并计算属性 / (新)对母网络进行稳健性模拟
  
  # (*** 新增: 对当前处理组的母网络进行稳健性模拟 ***)
  cat("  - 正在为该组的母网络计算稳健性...\n")
  N_SIMULATIONS_ROBUSTNESS <- 100 # 可以设置一个独立的模拟次数
  for (i in 1:N_SIMULATIONS_ROBUSTNESS) {
    robustness_value <- rand.remov.once(
      net_matrix = mother_adj_matrix,
      sp_ra = sp_ra
    )
    robustness_df <- data.frame(
      Proportion_removed = 0.5, # 固定移除比例
      Type = "Weight", # 固定类型
      group = current_treatment,
      reps = paste0("V", i), # 模拟重复ID
      values = robustness_value
    )
    all_robustness_results[[paste(current_treatment, i, sep="_")]] <- robustness_df
  }
  
  # (*** 原有逻辑: 提取子网并计算拓扑属性 ***)
  for (sample_id in colnames(otu_table_treatment)) {
    species_in_sample <- rownames(otu_table_treatment)[otu_table_treatment[, sample_id] > 0]
    nodes_for_subgraph <- intersect(species_in_sample, V(mother_network)$name)
    if (length(nodes_for_subgraph) > 1) {
      subnetwork_obj <- induced_subgraph(mother_network, vids = nodes_for_subgraph)
      if (inherits(subnetwork_obj, "igraph") && gorder(subnetwork_obj) > 0) {
        
        properties_list <- net_properties.4(subnetwork_obj)
        properties <- as.data.frame(t(properties_list))
        properties$vulnerability <- calculate_vulnerability(subnetwork_obj)
        properties <- as.data.frame(lapply(properties, function(x) as.numeric(as.character(x))))
        properties$SampleID <- sample_id
        all_properties_results[[paste(current_treatment, sample_id, sep="_")]] <- properties
      }
    }
  }
}

# --- 3. 格式化并保存 network_properties.csv ---
cat("\n--- 正在格式化并保存网络属性文件 ---\n")
if (length(all_properties_results) > 0) {
  combined_properties <- bind_rows(all_properties_results)
  properties_long <- combined_properties %>%
    rename(Samples = SampleID) %>%
    pivot_longer(cols = -Samples, names_to = "properties", values_to = "value") %>%
    mutate(properties = case_when(
      properties == "num.vertices.n." ~ "Nodes",
      properties == "num.edges.L." ~ "Links",
      properties == "average.degree.Average.K." ~ "Average_K",
      properties == "average.path.length" ~ "Average_path_length",
      properties == "no.clusters" ~ "Connected_components",
      properties == "RM.relative.modularity." ~ "RM",
      properties == "vulnerability" ~ "Vulnerability",
      TRUE ~ properties 
    )) %>%
    filter(properties %in% c("Links", "Nodes", "Average_K", "Average_path_length",
                             "Connected_components", "RM", "Vulnerability")) %>%
    select(Samples, properties, value)
  
  write.csv(properties_long, file.path("output", "network_properties.csv"), row.names = FALSE, na = "")
  cat("文件 'network_properties.csv' 已成功保存！\n")
} else {
  cat("警告: 未能计算任何网络属性，无法生成文件。\n")
}

# --- 4. 格式化并保存 robustness.csv ---
cat("\n--- 正在格式化并保存稳健性文件 ---\n")
if (length(all_robustness_results) > 0) {
  combined_robustness <- bind_rows(all_robustness_results) %>%
    select(Proportion_removed, Type, group, reps, values) # 调整列顺序以匹配示例
  
  write.csv(combined_robustness, file.path("output", "robustness.csv"), row.names = FALSE)
  cat("文件 'robustness.csv' 已成功保存！\n")
} else {
  cat("警告: 未能计算任何网络稳健性，无法生成文件。\n")
}

cat("\n所有任务完成！\n")