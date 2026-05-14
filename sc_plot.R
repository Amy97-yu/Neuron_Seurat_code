scplot = function(sce,features,idents) {

  library(Seurat)
  library(patchwork)
  library(tidyverse)
  library(ggsci)
  if (!dir.exists("./plot")) {
    dir.create("./plot")
  }
  #sce = sce
  #features <- read.table("./Genelist.txt")$V1
  #source("./hsa2mmu.R")
  #features <- hsa2mmu(features)
  #write.table(features,"gl_mmu.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)

  #features <- read.table("./gl_mmu.txt")$V1
  sce <- readRDS("./25.9.21_Seurat_Object_Annotated.rds")
  DefaultAssay(sce) <- "RNA"
  #Idents(sce) <- "seurat_clusters"
  sce[["RNA"]] <- JoinLayers(sce[["RNA"]])
  CellMarkers <- FindAllMarkers(sce,logfc.threshold = 0.25,only.pos = T)
  CellMarkers <- transform(CellMarkers,pct_FD = (pct.1 - pct.2))
  CellMarkers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 0.5) %>%
    slice_head(n = 10) %>% ungroup() -> top10
  DoHeatmap(sce, features = top10$gene,assay = "RNA",size = 5)
 source("E:/rnaseq/Seurat_PIP/Seurat_PIP/feature_plot_density.R") # 导入函数
  Genelist <- read.table("./NEU.txt")$V1

  for (i in Genelist) {

    p <- plot_density(obj = sce,
                      marker= i,
                      dim = "UMAP", size = 2,reduction = "umap",ncol = 3)
    p
    ggsave(paste("./plot/",i,"_plot_desity.pdf",sep = ""),height = 10,width = 10)
    ggsave(paste("./plot/",i,"_plot_desity.png",sep = ""),height = 10,width = 10)
  }

  write.csv(top5,"./top5_Cellmarker.csv")
  write.csv(CellMarkers,"./All_Cellmarker.csv")

  source("E:/rnaseq/Seurat_PIP/Seurat_PIP/sc_color_more.R")
  col <- scCOL(length(unique(Idents(sce)))) #画版
  group_col = scCOL(length(unique(sce$group)))
  ################################################################################

  # Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
  #DefaultAssay(sce) <- "RNA"
  #RidgePlot(sce, features = features[1], ncol = 1,cols = col)
  #ggsave("./plot/ridges_plot.pdf",width = 10,height = 10)
  #ggsave("./plot/ridges_plot.png",width = 10,height = 10)

  features <- read.table("./NEU.txt")$V1
  sce_feature <- sce[features,]

  features <- rownames(sce_feature)

  for (i in features) {

    p = VlnPlot(sce_feature, features = i,cols = col,pt.size = 1)

    p <- VlnPlot(sce_feature, features = i,cols = col,pt.size = 1)
    p
    ggsave(paste("./plot/",i,"_vlnplot.pdf",sep = ""),height = 5,width = 11)
    ggsave(paste("./plot/",i,"_vlnplot.png",sep = ""),height = 5,width = 11)
  }

  # Feature plot - visualize feature expression in low-dimensional space
  FeaturePlot(sce, features = features,cols = c("gray","red"),pt.size = 1,label = F)
  FeaturePlot(sce, 
              features = features,
              cols = c("lightgrey", "#E64B35"),  # 使用更专业的颜色搭配
              pt.size = 0.8,                     # 稍小的点尺寸避免过度重叠
              order = TRUE,                      # 表达高的点显示在上层
              blend = FALSE,
              combine = ifelse(length(features) > 1, TRUE, FALSE),  # 自动判断是否合并多图
              label = FALSE,
              label.size = 4,
              repel = TRUE) +                    # 避免标签重叠
    theme_classic() +                            # 使用更简洁的主题
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),  # 标题居中加粗
          legend.position = "right",            # 图例位置
          legend.title = element_text(size = 10),  # 图例标题
          legend.text = element_text(size = 8),  # 图例文字
          axis.line = element_line(color = "black"),  # 坐标轴线
          axis.text = element_text(color = "black", size = 10)) +  # 坐标轴文字
    guides(color = guide_colorbar(barwidth = 0.8, barheight = 4))  # 调整颜色条尺寸
  # Marker
  Idents(sce) <- "celltype"

  DotPlot(sce, features = features) + RotatedAxis()+ ggmin::theme_powerpoint()

  DotPlot(sce, features = features) + coord_flip()+
    theme_bw()+
    theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
    labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
    scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))

  DotPlot(sce, features = features) + 
    coord_flip() +
    theme_bw(base_size = 16) +  # 全局基础字体放大到16pt
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(
        size = 14,              # X轴文字大小
        hjust = 0.5,            # 关键修改：水平居中
        vjust = 0.5             # 垂直居中
      ),
      axis.text.y = element_text(size = 14),  # Y轴文字大小
      legend.text = element_text(size = 13),  # 图例文字大小
      legend.title = element_text(size = 14), # 图例标题大小
      plot.title = element_text(size = 16)    # 标题大小
    ) +
    labs(x = NULL, y = NULL) +
    guides(size = guide_legend(order = 8)) +
    scale_color_gradientn(
      values = seq(0, 1, 0.2),
      colours = c('#330066', '#336699', '#66CC66', '#FFCC33')
    )
  
  
  
  

  # Single cell heatmap of feature expression
  DefaultAssay(sce) <- "RNA"
  #Idents(sce) <- "celltype"
  DoHeatmap(subset(sce, downsample = 100), features = features, size = 3,slot = "data")
  DoHeatmap(subset(sce, downsample = 100), features = features, size = 3)

  DoHeatmap(subset(sce, downsample = 100), features = top5$gene, size = 3)

  pdf("./Pheatmap_top5.pdf",width = 10,height = 10,onefile = FALSE)
  DoHeatmap(subset(sce, downsample = 100), features = top5$gene, size = 3)
  dev.off()


  # Dimplot
  DimPlot(sce,label = T,repel = T,pt.size = 3,cols = col,reduction = "umap") +
            ggmin::theme_powerpoint()

  # 主题
  # remotes::install_github('sjessa/ggmin')
  DimPlot(sce,label = T,repel = T,pt.size = 3,cols = col) +
    ggmin::theme_powerpoint()


  sce_list <- SplitObject(sce,split.by = "group")

  p1 = DimPlot(sce_list[[1]],label = T,repel = T,pt.size = 3,cols = col) +
    ggmin::theme_powerpoint() + ggtitle(label = "D0") + theme(plot.title = element_text(hjust = 0.5))

  p2 = DimPlot(sce_list[[2]],label = T,repel = T,pt.size = 3,cols = col) +
    ggmin::theme_powerpoint() + ggtitle(label = "D1") + theme(plot.title = element_text(hjust = 0.5))

  p3 = DimPlot(sce_list[[3]],label = T,repel = T,pt.size = 3,cols = col) +
    ggmin::theme_powerpoint() + ggtitle(label = "D3") + theme(plot.title = element_text(hjust = 0.5))
 
  p4 = DimPlot(sce_list[[4]],label = T,repel = T,pt.size = 3,cols = col) +
    ggmin::theme_powerpoint() + ggtitle(label = "D7") + theme(plot.title = element_text(hjust = 0.5)) 
  

  (p1 + p2) / (p3 + p4) + plot_layout(guides = 'collect') + plot_layout(widths = c(2,2))
  return(sce)
}

