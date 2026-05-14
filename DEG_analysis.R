sce <- readRDS("../../../12_Seurat.rds")
dir.create("../../../13_IntersetGene")
setwd("../../../13_IntersetGene/")
# DEGs MT  sigGene
DEGs = features
length(DEGs)
#unique(DEGs)
MT <- read.table("../mt/MT.txt")$V1
sigGenes <- readRDS("../10_monocle/monocle2_v2.1/monocle2_v2.1/Result_0303/sig.gene.rds")


tmp = intersect(sigGenes,DEGs)
write(tmp,"./tmp_sigGenes&DEGs.txt")


#install.packages("ggvenn")
library(ggvenn)

gl_list = list(DEGs = DEGs,MT = MT,sigGenes = sigGenes)
ggvenn::ggvenn(data = gl_list,columns = c("DEGs","MT","sigGenes"))

sl_genes = intersect(intersect(DEGs,MT),sigGenes)
length(sl_genes)
write(sl_genes,"sl_genes.txt")

DEG_anaysis <- function(sce,idents_celltype) {
  library(Seurat)
  library(tidyverse)
  library(patchwork)

  if (!dir.exists("./DEGs_analysis")) {
    dir.create("./DEGs_analysis")
  }
  setwd("./DEGs_analysis")
  sce <- readRDS("../2024.12.14_sce_anno_FINAL.rds")
  Idents(sce) <- "group"
  # Idents(sce) <- "time"      # DEG  比较的是时间之间的差异 ===  bulkRNAseq  sce all celltype
  DEGs <- FindAllMarkers(sce,logfc.threshold = 1,only.pos = FALSE,min.pct = 0.3) # 差异基因上下调都保留

  DEGs  <- DEGs[which(DEGs$p_val_adj <= 0.01),]

  features <- DEGs$gene
  length(features)


  DEGs %>% group_by(cluster) %>% filter(avg_log2FC > 1) %>%
    slice_head(n = 30) %>% ungroup() -> top10_DEGs
  source("E:/rnaseq/Seurat_PIP/Seurat_PIP/sc_plot.R")
  features <- top10_DEGs$gene
  sce <- scplot(sce = sce,features = features,idents = "group")
  ################################################################################
  #idents_celltype <- "seurat_clusters"
  Idents(sce) <- "celltype"
  cellnames <- unique(sce@active.ident)

  for (i in cellnames) {
    #i = cellnames[1]
    Result_dir = paste("Result",i,sep = "-")

    if (!dir.exists(Result_dir)) {
      dir.create(Result_dir)
    }

    setwd(Result_dir)
    sce$cellname <- Idents(sce)
    sce_sub <- subset(sce,cellname %in% i)

    Idents(sce_sub) <- "group"
    # Idents(sce) <- "time"      # DEG  比较的是时间之间的差异 ===  bulkRNAseq  sce all celltype
    DEGs <- FindAllMarkers(sce_sub,logfc.threshold = 0.25,only.pos = TRUE) # 差异基因上下调都保留
    DEGs %>% group_by(cluster) %>% filter(avg_log2FC > 1) %>%
      slice_head(n = 12) %>% ungroup() -> top12
    DoHeatmap(subset(sce_sub, downsample = 100), features = top12$gene, size = 3,slot = "data")
    write.csv(DEGs,"../12.15.markers_DEGs_time.csv")
    source("E:/rnaseq/Seurat_PIP/Seurat_PIP/sc_plot.R")
    features <- top10_DEGs$gene
    sce_sub <- scplot(sce = sce_sub,features = features,idents = "group")
    setwd("../")
  }
  ################################################################################
  return(sce)
}
