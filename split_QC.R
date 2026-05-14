##########分开做质控############################################################
split_QC <- FALSE #TRUE
if (split_QC) {
  sce_list <- SplitObject(sce,split.by = "sample")

  sce_list$D0$percent.mt <- Seurat::PercentageFeatureSet(sce_list$D0, pattern = "^mt")
  # 小鼠细胞中红细胞基因主要有哪些？
  # UMI Gene
  VlnPlot(sce_list$D0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0,group.by = "group")
  sce_list$D0 <- subset(sce_list$D0, subset = nFeature_RNA > 2500 & nFeature_RNA < 7000 & percent.mt < 20)
  VlnPlot(sce_list$D0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0,group.by = "group")
}
sce_list$D1$percent.mt <- Seurat::PercentageFeatureSet(sce_list$D1, pattern = "^mt")
# 小鼠细胞中红细胞基因主要有哪些？
# UMI Gene
VlnPlot(sce_list$D, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0,group.by = "group")
sce_list$D0 <- subset(sce_list$D0, subset = nFeature_RNA > 2500 & nFeature_RNA < 7000 & percent.mt < 20)
VlnPlot(sce_list$D0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0,group.by = "group")
}
################################################################################
