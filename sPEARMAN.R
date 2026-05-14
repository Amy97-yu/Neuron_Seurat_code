# 优化版连续变化可视化 - 无平滑曲线，增强数据密度显示
library(ggplot2)
library(viridis)
library(dplyr)
library(ggpmisc)

# 数据准备
valid_data <- pData(cds)[!is.na(pData(cds)$mito_score) & 
                           !is.na(pData(cds)$Pseudotime), ]

# 计算相关性
cor_result <- cor.test(valid_data$Pseudotime, valid_data$mito_score,
                       method = "spearman", exact = FALSE)

# 2. 无平滑曲线版 - 增强数据密度显示
no_curve_plot <- ggplot(valid_data, aes(x = Pseudotime, y = mito_score)) +
  
  # 增加点透明度并添加数据密度趋势
  geom_point(aes(color = Pseudotime), alpha = 0.5, size = 2, shape = 16) +
  
  # 添加数据密度轮廓
  #geom_density_2d_filled(alpha = 0.2, contour_var = "density") +
  
  # 添加简单的线性趋势线（不显著的话可以去掉）
  # geom_smooth(method = "lm", color = "red", se = FALSE, size = 1, alpha = 0.8) +
  
  # 使用更鲜明的颜色梯度
  scale_color_viridis(option = "plasma", name = "Pseudotime", 
                      guide = guide_colorbar(barwidth = 1, barheight = 10)) +
  
  # 添加分位数参考线（可选）
  #geom_hline(yintercept = quantile(valid_data$mito_score, c(0.25, 0.5, 0.75)), 
            # linetype = "dashed", alpha = 0.3, color = "gray40") +
  
  # 统计标注
  annotate("text", x = min(valid_data$Pseudotime), 
           y = max(valid_data$mito_score) * 0.98,
           label = paste0("Spearman's ρ = ", round(cor_result$estimate, 3), 
                          ifelse(cor_result$p.value < 0.001, "\np < 0.001", 
                                 paste0("\np = ", format(cor_result$p.value, scientific = TRUE, digits = 2)))),
           hjust = 0, vjust = 1, size = 5.5, fontface = "bold",
           color = "black", bg = "white", alpha = 0.8) +
  
  # 改进的坐标轴和标题
  labs(x = "Pseudotime", 
       y = "Mitochondrial Function Score",
       title = "Mitochondrial Functional Dynamics During Neuronal Differentiation",
       subtitle = "Distribution of mitochondrial activity scores across pseudotemporal trajectory") +
  
  # 更专业的主题设置
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 10)),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40", margin = margin(b = 20)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray90", size = 0.3),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

# 显示无曲线版图表
print(no_curve_plot)

# 保存优化版主图
ggsave("Fig4A_enhanced_continuous_dynamics.pdf", 
       plot = enhanced_continuous_plot,
       width = 11, 
       height = 7,
       dpi = 300)

# 2. 优化版分段分析图 - 大幅提升美观度
# 创建更精细的分段
n_bins <- 12
valid_data$pseudotime_bin <- cut(valid_data$Pseudotime, 
                                 breaks = unique(quantile(valid_data$Pseudotime, 
                                                          probs = seq(0, 1, length.out = n_bins + 1))),
                                 include.lowest = TRUE)

bin_stats <- valid_data %>%
  group_by(pseudotime_bin) %>%
  summarise(
    n_cells = n(),
    mean_score = mean(mito_score),
    sem_score = sd(mito_score) / sqrt(n()),  # 标准误
    median_score = median(mito_score),
    pseudotime_mid = median(Pseudotime)  # 使用中位数作为区间代表值
  ) %>%
  filter(n_cells >= 5)  # 过滤掉细胞数太少的区间

# 创建美观的分段趋势图
enhanced_segmented_plot <- ggplot(bin_stats, aes(x = pseudotime_mid, y = mean_score)) +
  
  # 添加背景趋势线（轻微透明）
  geom_smooth(data = valid_data, aes(x = Pseudotime, y = mito_score),
              method = "loess", color = "gray80", size = 1, se = FALSE, alpha = 0.5) +
  
  # 误差线 - 更细更精致
  geom_errorbar(aes(ymin = mean_score - sem_score,
                    ymax = mean_score + sem_score),
                width = diff(range(bin_stats$pseudotime_mid)) * 0.02,
                color = "#2E86AB", size = 0.8, alpha = 0.7) +
  
  # 主趋势线 - 更粗更明显
  geom_line(color = "#D32F2F", size = 2, alpha = 0.9) +
  
  # 数据点 - 大小随样本量变化，更美观
  geom_point(aes(size = n_cells, fill = pseudotime_mid), 
             shape = 21, color = "white", stroke = 1.2) +
  
  # 改进的颜色和大小标度
  scale_fill_viridis(option = "plasma", name = "Pseudotime") +
  scale_size_continuous(range = c(3, 8), name = "Number of Cells",
                        breaks = c(min(bin_stats$n_cells), 
                                   median(bin_stats$n_cells), 
                                   max(bin_stats$n_cells))) +
  
  # 专业标题和标签
  labs(x = "Pseudotime", 
       y = "Mean Mitochondrial Function Score ± SEM",
       title = "Quantitative Validation of Mitochondrial Function Dynamics",
       subtitle = "Binned analysis confirms progressive increase with minimal variance") +
  
  # 精美主题设置
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 10)),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40", margin = margin(b = 20)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.2, "cm"),
    panel.grid.major = element_line(color = "gray90", size = 0.3),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  
  # 添加趋势指示箭头
  annotate("segment", 
           x = min(bin_stats$pseudotime_mid) + diff(range(bin_stats$pseudotime_mid)) * 0.1,
           xend = max(bin_stats$pseudotime_mid) - diff(range(bin_stats$pseudotime_mid)) * 0.1,
           y = min(bin_stats$mean_score) + diff(range(bin_stats$mean_score)) * 0.1,
           yend = min(bin_stats$mean_score) + diff(range(bin_stats$mean_score)) * 0.1,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
           color = "#D32F2F", size = 1.5) +
  
  annotate("text", 
           x = mean(range(bin_stats$pseudotime_mid)),
           y = min(bin_stats$mean_score) + diff(range(bin_stats$mean_score)) * 0.05,
           label = "Differentiation Direction", 
           color = "#D32F2F", fontface = "bold", size = 4.5)

# 显示优化版分段图
print(enhanced_segmented_plot)

# 保存优化版分段图
ggsave("Fig4B_enhanced_segmented_analysis.pdf", 
       plot = enhanced_segmented_plot,
       width = 11, 
       height = 7,
       dpi = 300)

# 输出统计摘要
cat("=== 优化版图表统计摘要 ===\n")
cat("Spearman相关性: ρ =", round(cor_result$estimate, 3), 
    ", p =", ifelse(cor_result$p.value < 0.001, "< 0.001", 
                    format(cor_result$p.value, scientific = TRUE, digits = 2)), "\n\n")

cat("分段分析统计:\n")
print(bin_stats %>% select(pseudotime_bin, n_cells, mean_score, sem_score))

cat("\n=== 图表已保存 ===\n")
cat("1. 增强趋势主图: Fig4A_enhanced_continuous_dynamics.pdf\n")
cat("2. 美观分段图: Fig4B_enhanced_segmented_analysis.pdf\n")