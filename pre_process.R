setwd("./14_sc_sop/")

pre_process = function(data_dir,sample_info) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(tidyverse)
    library(patchwork)
  })

  data_dir <- "../../01_RawData/data/"
  sample_info <- data.frame(sample = dir(data_dir),group = c("D0","D1","D3","D7"))

  D_list <- list()
  for (i in 1:length(sample_info$sample)) {
    D_data <- Read10X(paste(data_dir,sample_info$sample[i],sep = ""))
    D_sce <- CreateSeuratObject(D_data,min.cells = 3,min.features = 200)
    D_sce@meta.data$group <- sample_info$group[i]
    D_sce@meta.data$sample <- sample_info$sample[i]
    D_list = c(D_list,list(D_sce))
  }
  sce <- merge(D_list[[1]],D_list[2:length(D_list)],add.cell.ids = sample_info$sample)
  sce_list <- SplitObject(sce,split.by = "group")
  D_list <- sce_list

  D_list <- lapply(X = D_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })

  ifnb.list <- D_list
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = ifnb.list)
  ################################################################################
  immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features) # 锚点
  immune.combined <- IntegrateData(anchorset = immune.anchors)

  # 线粒体含量
  DefaultAssay(immune.combined) <- "RNA"
  Idents(immune.combined) <- "group"
  immune.combined$percent.mt <- Seurat::PercentageFeatureSet(immune.combined, pattern = "^mt") # MT
  # 小鼠细胞中红细胞基因主要有哪些？
  # UMI Gene
  VlnPlot(immune.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0,group.by = "group")

  immune.combined <- subset(immune.combined, subset = nFeature_RNA > 2500 & nFeature_RNA < 7000 & percent.mt < 10) # 400 - 7000
  VlnPlot(immune.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0,group.by = "group")
  saveRDS(immune.combined,"./All_sce_filter_FINAL.rds")
  return(immune.combined) # 返回的对象sce对象
  ################################################################################
}
# 首先检查各样本的基因数量和基因名是否一致
gene_counts <- sapply(D_list, function(x) nrow(x))
gene_names <- lapply(D_list, function(x) rownames(x))

cat("各样本基因数量:\n")
print(gene_counts)

# 检查是否有基因数量不一致的样本
if (length(unique(gene_counts)) > 1) {
  cat("警告: 样本间基因数量不一致!\n")
  
  # 找出所有样本共有的基因
  common_genes <- Reduce(intersect, gene_names)
  cat("共有基因数量:", length(common_genes), "\n")
  
  # 只保留共有基因
  D_list_filtered <- lapply(D_list, function(x) {
    x <- x[common_genes, ]
    return(x)
  })
  
  # 使用过滤后的基因列表重新进行处理
  D_list <- lapply(D_list_filtered, function(x) {
    # 数据标准化
    x <- NormalizeData(
      object = x,
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    )
    
    # 寻找高变基因
    x <- FindVariableFeatures(
      object = x,
      selection.method = "vst",
      nfeatures = min(2000, nrow(x)),  # 确保不超过基因总数
      verbose = FALSE
    )
    
    return(x)  # 必须返回处理后的对象
  })
  
} else {
  # 如果基因数量一致，尝试单独处理每个样本
  D_list <- lapply(D_list, function(x) {
    tryCatch({
      # 数据标准化
      x <- NormalizeData(
        object = x,
        normalization.method = "LogNormalize",
        scale.factor = 10000,
        verbose = FALSE
      )
      
      # 寻找高变基因
      x <- FindVariableFeatures(
        object = x,
        selection.method = "vst",
        nfeatures = min(2000, nrow(x)),
        verbose = FALSE
      )
      
      return(x)
    }, error = function(e) {
      cat("处理样本时出错:", e$message, "\n")
      return(x)  # 返回未处理的样本
    })
  })
}

# 检查处理结果
cat("处理完成。各样本信息:\n")
for (i in 1:length(D_list)) {
  cat("样本", i, ":")
  cat(" 细胞数:", ncol(D_list[[i]]))
  cat(" 基因数:", nrow(D_list[[i]]))
  cat(" 高变基因数:", length(VariableFeatures(D_list[[i]])), "\n")
}

# 如果仍有问题，可以尝试逐个样本处理
# 这样可以更清楚地看到哪个样本出了问题
for (i in 1:length(D_list)) {
  cat("正在处理样本", i, "...\n")
  tryCatch({
    D_list[[i]] <- NormalizeData(D_list[[i]], verbose = FALSE)
    D_list[[i]] <- FindVariableFeatures(
      D_list[[i]], 
      selection.method = "vst", 
      nfeatures = 2000,
      verbose = FALSE
    )
    cat("样本", i, "处理成功\n")
  }, error = function(e) {
    cat("样本", i, "处理失败:", e$message, "\n")
  })
}
library(ggplot2)
library(Seurat)
library(ggsci)  # 用于Nature风格配色

# 设置Nature风格主题
nature_theme <- theme(
  plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
  axis.title = element_text(size = 14, face = "bold"),
  axis.text = element_text(size = 12, color = "black"),
  axis.text.x = element_text(angle = 45, hjust = 1),
  legend.title = element_text(size = 12, face = "bold"),
  legend.text = element_text(size = 10),
  panel.background = element_blank(),
  panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
  strip.background = element_rect(fill = "grey90", color = "black"),
  strip.text = element_text(size = 12, face = "bold")
)

# 过滤前的质控图
p1 <- VlnPlot(immune.combined, 
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
              ncol = 3, 
              pt.size = 0, 
              group.by = "group") +
  scale_fill_npg() +  # Nature风格配色
  nature_theme +
  ggtitle("Pre-filtering Quality Metrics") +
  theme(plot.title = element_text(hjust = 0.5))

print(p1)

# 应用过滤条件
immune.combined <- subset(immune.combined, 
                          subset = nFeature_RNA > 2500 & 
                            nFeature_RNA < 7000 & 
                            percent.mt < 10)

# 过滤后的质控图
p2 <- VlnPlot(immune.combined, 
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
              ncol = 3, 
              pt.size = 0, 
              group.by = "group") +
  scale_fill_npg() +  # Nature风格配色
  nature_theme +
  ggtitle("Post-filtering Quality Metrics") +
  theme(plot.title = element_text(hjust = 0.5))

print(p2)

# 组合两个图表以便比较
library(patchwork)
combined_plot <- p1 / p2 + 
  plot_annotation(title = "Quality Control Metrics Before and After Filtering",
                  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))

print(combined_plot)

# 保存高质量图片
ggsave("qc_metrics_nature_style.png", 
       plot = combined_plot, 
       width = 12, 
       height = 10, 
       dpi = 300)

# 可选：添加统计信息到图表
# 计算各组的统计信息
qc_stats <- immune.combined@meta.data %>%
  group_by(group) %>%
  summarise(
    mean_nFeature = mean(nFeature_RNA),
    mean_nCount = mean(nCount_RNA),
    mean_mt = mean(percent.mt),
    cells_remaining = n()
  )

print(qc_stats)

# 将统计信息添加到图表中（可选）
p2_with_stats <- p2 +
  geom_text(data = qc_stats, 
            aes(x = group, y = max(immune.combined$nFeature_RNA) * 1.1, 
                label = paste0("n=", cells_remaining)), 
            size = 4, fontface = "bold")

print(p2_with_stats)
