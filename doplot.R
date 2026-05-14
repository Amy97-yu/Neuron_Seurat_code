library(SeuratData)
library(patchwork)
library(dbplyr)
library(tidyverse)
library(Seurat)
sce <- readRDS("../25.9.20 chongxinfenxi/25.9.21_Seurat_Object_Annotated.rds")
DimPlot(sce, reduction = "umap",label = T)
# 按类别排序的标记物列
markers_to_check <- c(
  # 干细胞/祖细胞标记
  "Nes", "Gfap", "rna_Sox2", "Pax6",
  # 神经元标记
  "Ascl1", "Nrxn1", "rna_Map2",
  # 星形胶质细胞标记
  "Gja1", "Aqp4", "Agt",
  # 少突胶质细胞谱系标记
  "Cspg4", "Olig2", "Vcan"
)

# 强制指定细胞簇的顺序 - 请根据您的实际情况调整这些簇名
# 这里假设您的细胞簇按照神经干细胞→神经元→星形胶质细胞→少突胶质细胞的顺序排列
desired_cluster_order <- c("Neural Stem Cells", "NPCs","Neurons", "Astrocytes", "OPCs")
# 如果您的簇名不同，请替换为实际的簇名

# 设置细胞簇顺序
Idents(sce) <- factor(Idents(sce), levels = desired_cluster_order)

# 生成更醒目的点图
p <- DotPlot(sce, features = markers_to_check, 
             dot.scale = 8,  # 增大点的大小
             cols = c("#2E86AB", "#A23B72")) +  # 使用更鲜明的颜色
  RotatedAxis() +
  scale_color_gradientn(colors = viridis::viridis(20)) +  # 使用更丰富的颜色渐变
  labs(x = "Gene Markers", 
       y = "Cell Types", 
       title = "Marker Expression Profile") +
  theme_minimal(base_size = 14) +  # 增大基础字体大小
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    panel.grid.major = element_line(color = "grey80"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  )

print(p) 
