rm(list = ls())
pro_dir <- "./03_analysis/14_sc_sop/"
################################################################################
source("./pre_process.R") #运行一个脚本
data_dir <- "D:/Work/03_NyProject/03_教学/scRNA_yuan/03_analysis/03_Conditon/data/"
sample_info <- data.frame(sample = dir(data_dir),group = c("D0","D1","D3","D7"))
sce <- pre_process(data_dir = data_dir,sample_info = sample_info)
################################################################################
source("./Seurat_pip.R")
sce <- sce
res <- 2
sce <- Seurat_pip(sce = sce,res = res)
################################################################################
# 输出表格
# 绘图
source("./sc_plot.R")
features <- read.table("./gl_mmu.txt")$V1
sce <- scplot(sce = sce,features = features,idents = "seurat_clusters")
################################################################################
#setwd("../../")
source("./DEG_analysis.R")
sce <- DEG_anaysis(sce = sce,idents_celltype = "seurat_clusters")
################################################################################
