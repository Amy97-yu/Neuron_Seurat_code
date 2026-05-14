# 定义氧化磷酸化核心基因集
oxphos_genes <- c("Abcb7", "Acaa1", "Acaa2", "Acdam", "Acadsb", "Acadvl", "Acat1", "Aco2", "Afg3l2", "Aifm1", "Alas1", "Aldh6a1", "Atp1b1", "Atp6ap1", "Atp6v0b", "Atp6v0c", "Atp6v0e1", "Atp6v1c1", "Atp6v1d", "Atp6v1e1", "Atp6v1f", "Atp6v1g1", "Atp6v1h", "Bax", "Bckdha", "Bdh2", "Casp7", "Cox10", "Cox11", "Cox15", "Cox17", "Cox4i1", "Cox5a", "Cox6a1", "Cox6b1", "Cox7a2", "Cox7a2l", "Cox7b", "Cox8a", "Cpt1a", "Cs", "Cyc1", "Decr1", "Dlat", "Dld", "Dlst", "Ech1", "Echs1", "Eci1", "Etfa", "Etfb", "Etfdh", "Fdx1", "Fh", "Fxn", "Glud1", "Got2", "Gpi", "Gpx4", "Grpel1", "Hadha", "Hadhb", "Hccs", "Hsd17b10", "Hspa9", "Htra2", "Idh1", "Idh2", "Idh3a", "Idh3b", "Idh3g", "Immt", "Iscu", "Ldha", "Ldhb", "Lrpprc", "Maob", "Mdh1", "Mdh2", "Mfn2", "Mgst3", "Mpc1", "Mrpl11", "Mrpl15", "Mrpl34", "Mrpl35", "Mrps11", "Mrps12", "Mrps15", "Mrps22", "Mrps30", "Mtrf1", "Mtrr", "Mtx2", "Ndufa1", "Ndufa2", "Ndufa3", "Ndufa4", "Ndufa5", "Ndufa6", "Ndufa7", "Ndufa8", "Ndufa9", "Ndufab1", "Ndufb2", "Ndufb3", "Ndufb5", "Ndufb6", "Ndufb7", "Ndufb8", "Ndufc1", "Ndufc2", "Ndufs1", "Ndufs2", "Ndufs3", "Ndufs4", "Ndufs6", "Ndufs7", "Ndufs8", "Ndufv1", "Ndufv2", "Nqo2", "Oat", "Ogdh", "Opa1", "Oxa1l", "Pdha1", "Pdhb", "Pdhx", "Pdk4", "Pdp1", "Phb2", "Phyh", "Pmpca", "Polr2f", "Por", "Prdx3", "Retsat", "Rhot1", "Rhot2", "Sdha", "Sdhb", "Sdhc", "Sdhd", "Slc25a11", "Slc25a12", "Slc25a20", "Slc25a3", "Slc25a4", "Slc25a5", "Sucla2", "Suclg1", "Supv3l1", "Surf1", "Tcirg1", "Timm10", "Timm13", "Timm17a", "Timm50", "Timm8b", "Timm9", "Tomm22", "Tomm70a", "Uqcr10", "Uqcr11", "Uqcrb", "Uqcrc1", "Uqcrc2", "Uqcrfs1", "Uqcrhl", "Uqcrq", "Vdac1", "Vdac2", "Vdac3", "Ndufa10", "Cycs", "Cox5b", "Cox6c", "Cox7c", "Atp5a1", "Atp5b", "Atp5c1", "Atp5d", "Atp5e", "Atp5f1", "Atp5g1", "Atp5g2")

# 检查基因在数据中的存在情况
available_genes <- oxphos_genes[oxphos_genes %in% rownames(exprs(cds))]
print(paste("Found", length(available_genes), "OXPHOS genes (out of", length(oxphos_genes), ")"))

# 安装和加载必要包
if (!require("AUCell", quietly = TRUE)) {
  BiocManager::install("AUCell")
}
if (!require("viridis", quietly = TRUE)) {
  install.packages("viridis")
}
library(AUCell)
library(viridis)
library(ggplot2)
library(ggpubr)

# 提取表达矩阵
expr_matrix <- as.matrix(exprs(cds))
print(paste("Expression matrix dimensions:", dim(expr_matrix)[1], "genes ×", dim(expr_matrix)[2], "cells"))

# 计算AUCell评分
cells_rankings <- AUCell_buildRankings(expr_matrix, 
                                       nCores = 1,
                                       plotStats = FALSE,
                                       verbose = FALSE)

cells_AUC <- AUCell_calcAUC(list(OXPHOS = available_genes), 
                            cells_rankings, 
                            aucMaxRank = ceiling(0.1 * nrow(cells_rankings)))

oxphos_scores <- getAUC(cells_AUC)["OXPHOS", ]
pData(cds)$mito_score <- as.numeric(oxphos_scores)

print("AUCell scoring completed!")
print(summary(pData(cds)$mito_score))

# 统计检验
if (!"Pseudotime" %in% colnames(pData(cds))) {
  stop("Error: Pseudotime column not found. Please complete pseudotime analysis first.")
}

cor_result <- cor.test(pData(cds)$Pseudotime, pData(cds)$mito_score, 
                       method = "spearman", exact = FALSE)

print("=== Statistical Test Results ===")
print(paste("Spearman's rho:", round(cor_result$estimate, 4)))
print(paste("P-value:", signif(cor_result$p.value, 3)))
print(paste("Significance:", ifelse(cor_result$p.value < 0.05, "Significant", "Not significant")))

# 创建更美观的可视化
main_plot <- ggplot(pData(cds), aes(x = Pseudotime, y = mito_score)) +
  geom_point(aes(color = mito_score), alpha = 0.7, size = 1.5) +
  geom_smooth(method = "loess", color = "#2E86AB", fill = "#A23B72", alpha = 0.2, size = 1.2) +
  scale_color_viridis(option = "plasma", name = "Mitochondrial\nfunction score") +
  labs(title = "Mitochondrial Function During Neuronal Differentiation",
       subtitle = paste("Spearman's ρ =", round(cor_result$estimate, 3), 
                        ", p =", signif(cor_result$p.value, 3)),
       x = "Pseudotime", 
       y = "Mitochondrial function score") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.major = element_line(color = "gray90", size = 0.2),
    panel.grid.minor = element_blank()
  )

# 按细胞类型分组显示（如果存在细胞类型信息）
if ("CellType" %in% colnames(pData(cds))) {
  # 定义美观的细胞类型颜色
  celltype_colors <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", 
                       "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7")
  
  celltype_plot <- ggplot(pData(cds), aes(x = CellType, y = mito_score, fill = CellType)) +
    geom_violin(alpha = 0.8, trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "white", alpha = 0.8, outlier.shape = NA) +
    stat_compare_means(method = "kruskal.test", 
                       label = "p.format", 
                       label.y = max(pData(cds)$mito_score) * 1.1) +
    scale_fill_manual(values = celltype_colors) +
    labs(title = "Mitochondrial Function Across Cell Types", 
         x = "Cell Type", 
         y = "Mitochondrial function score") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 组合图形
  combined_plot <- ggarrange(main_plot, celltype_plot, 
                             ncol = 2, 
                             widths = c(2, 1),
                             labels = c("A", "B"))
  print(combined_plot)
  
  # 保存高质量图片
  ggsave("mitochondrial_function_analysis.tiff", 
         plot = combined_plot,
         width = 12, 
         height = 5, 
         dpi = 300)
  
} else {
  print(main_plot)
  # 保存高质量图片
  ggsave("mitochondrial_function_analysis.tiff", 
         plot = main_plot,
         width = 8, 
         height = 6, 
         dpi = 300)
}

# 添加相关性文本输出
cor_text <- paste("Mitochondrial function score showed", 
                  ifelse(cor_result$estimate > 0, "positive", "negative"),
                  "correlation with pseudotime (Spearman's ρ =",
                  round(cor_result$estimate, 3), ", p =", 
                  ifelse(cor_result$p.value < 0.001, "p < 0.001", 
                         paste("p =", round(cor_result$p.value, 3))), ")")

print(cor_text)
