
Seurat_pip = function(sce,res) {
  sce <- readRDS("./All_sce_filter_FINAL.rds")
  immune.combined <- sce
  DefaultAssay(immune.combined) <- "integrated"
  immune.combined <- ScaleData(immune.combined, verbose = FALSE)
  immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
  immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
  immune.combined <- RunTSNE(immune.combined, reduction = "pca", dims = 1:30)
  immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:15)
  immune.combined <- FindClusters(immune.combined,resolution = 2.7)
  # Visualization
  DimPlot(immune.combined, reduction = "tsne",pt.size = 2,label = T)
  DimPlot(immune.combined, reduction = "umap", split.by  = "group",pt.size = 1,ncol = 2)
  DimPlot(immune.combined, reduction = "umap",pt.size = 1,label = T)
  ggsave("./umap_allcells.pdf",width = 10,height = 7)
  saveRDS(immune.combined,"./2025.5.29 All_sce_raw.rds")
  return(immune.combined)
}

